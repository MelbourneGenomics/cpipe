/////////////////////////////////////////////////////////////////////////////////
//
// This file is part of Cpipe.
// 
// Cpipe is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, under version 3 of the License, subject
// to additional terms compatible with the GNU General Public License version 3,
// specified in the LICENSE file that is part of the Cpipe distribution.
// 
// Cpipe is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Cpipe.  If not, see <http://www.gnu.org/licenses/>.
// 
/////////////////////////////////////////////////////////////////////////////////

load "pipeline_helpers.groovy"

///////////////////////////////////////////////////////////////////
// stages
///////////////////////////////////////////////////////////////////

check_sample_info = {

    doc "Validate basic sample information is correct"

    stage_status("check_sample_info", "enter", "n/a");

    def missingSummary = []
    for(sample in all_samples) {

        // Check that FASTQ files begin with the sample name followed by underscore
        def files = sample_info[sample].files.fastq
        if(files.any { !file(it).name.startsWith(sample+"_")}) {
            files.each { println "FASTQ: $it | sample=$sample" }
            fail report('templates/invalid_input.html') to channel: cpipe_operator, 
                                                           subject: "FASTQ files for sample $sample have invalid file name format", 
                                                           message: "Files $files do not start with the sample name $sample" 
        }

        // Check that all the files specified for the sample exist
        def missingFiles = files.grep { !file(it).exists() }
        if(missingFiles) 
            missingSummary << """
                    The following files specified for sample $sample could not be found:\n\n${missingFiles*.center(120).join('\n')}

                    Please check that the files in your sample file really exist in the data directory.
                 """.stripIndent()

        // Check that file names contain the lane information
        def missingLanes = files.grep { !(it ==~ ".*_L[0-9]*_.*") }
        if(missingLanes) 
            missingSummary << """
                    The following files specified for sample $sample do not contain lane information:\n\n${missingLanes*.center(120).join('\n')}

                    FASTQ file names are required to contain lane identifiers such as L001, L1 or similar. 
                    Please check your input FASTQ and rename it if necessary.
            """

        // Check that file names contain the read number information
        def missingRP = files.grep { !(it ==~ ".*_R[0-9][_.].*fastq.gz\$") }
        if(missingRP) 
            missingSummary << """
                    The following files for sample $sample do not contain the read number in the expected format:\n\n${missingRP*.center(120).join('\n')}

                    FASTQ file names are required to contain the number of the read from the read pair (1 or 2) 
                    in the form '_R1_' or '_R1.'. Please check your input FASTQ and rename it if necessary.
            """
    }

    if(missingSummary) {
        fail missingSummary.join("\n" + ("-" * 120) + "\n")
    }

    stage_status("check_sample_info", "exit", "n/a");
}

check_tools = {
    doc """
        Checks for presence of optional tools and sets appropriate pipeline variables
        to enable or disable corresponding pipeline features
        """

    stage_status("check_tools", "enter", "n/a");

    var UPDATE_VARIANT_DB : VARIANT_DB,
        ANNOTATION_VARIANT_DB : VARIANT_DB

    produce("revision.txt") {
        exec """
            git describe --always > $output.txt || true
        """
    }

    if(file(GROOVY_NGS).name in ["1.0.1","1.0"])
        fail "This version of Cpipe requires GROOVY_NGS >= 1.0.2. Please edit config.groovy to set the latest version of tools/groovy-ngs-utils"

    branch.UPDATE_VARIANT_DB = UPDATE_VARIANT_DB
    branch.ANNOTATION_VARIANT_DB = ANNOTATION_VARIANT_DB

    stage_status("check_tools", "exit", "n/a");
}

update_gene_lists = {
    doc "find additionally specified genes and add new ones that aren't on the incidentalome, to the gene lists"

    // builds additional genes from sample metadata file, then adds any new ones to the flagship
    // creates files: ../design/cohort.add.genes.txt, cohort.addonce.sample.genes.txt, cohort.notfound.genes.txt
    produce ('update_gene_lists.log') {
        from(sample_metadata_file) {
            exec """
                mkdir -p "../design"

                # This updates the local design's gene list with any genes in the prioritised genes that were missing from it
                python $SCRIPTS/find_new_genes.py --reference "$BASE/designs/genelists/exons.bed" --exclude "$BASE/designs/genelists/incidentalome.genes.txt" --target ../design < $sample_metadata_file

                # This script transfers this updated gene list to the global design, which is shared amongst all batches using it
                # This could be disabled because this might cause unpredictable behaviour
                python $SCRIPTS/update_gene_lists.py --source ../design --target "$BASE/designs" --log "$BASE/designs/genelists/changes.genes.log"

                touch update_gene_lists.log
            """
        }
    }
}

create_combined_target = {

    // Construct the region for variant calling from 
    //
    //   a) all of the disease cohort BED files
    //   b) the EXOME target regions
    //   c) any additional genes being analyzed
    //
    // This way we avoid calling variants over the entire genome, but still
    // include everything of interest
    String diseaseGeneLists = ANALYSIS_PROFILES.collect { "$BASE/designs/${it}/${it}.genes.txt" }.join(" ")
    String diseaseBedFiles = ANALYSIS_PROFILES.collect { "$BASE/designs/${it}/${it}.bed" }.join(" ")

    output.dir = "../design"

    produce("combined_target_regions.bed") {
        exec """
            python $SCRIPTS/combine_target_regions.py --genefiles $diseaseGeneLists --genefiles_required ../design/*.addonce.*.genes.txt --bedfiles $diseaseBedFiles $EXOME_TARGET --exons $BASE/designs/genelists/exons.bed |
            cut -f 1,2,3 | 
            $BEDTOOLS/bin/bedtools sort | 
            $BEDTOOLS/bin/bedtools merge > $output.bed
        """
    }

    branch.COMBINED_TARGET = output.bed
    exec """
        echo "create_combined_target: combined target is $COMBINED_TARGET"
    """
}

create_synonymous_target = {
    doc "find regions that allow synonymous variants"

    output.dir = "../design"
    produce( "combined_synonymous_regions.bed" ) {
        def safe_tmp_dir = [TMPDIR, UUID.randomUUID().toString()].join( File.separator )

        exec """
            mkdir -p "$safe_tmp_dir"

            $BEDTOOLS/bin/bedtools slop -i $input.bed -g $HG19_CHROM_INFO -b $ALLOW_SYNONYMOUS_INTRON > "$safe_tmp_dir/intron.bed"

            $BEDTOOLS/bin/bedtools slop -i $input.bed -g $HG19_CHROM_INFO -b -$ALLOW_SYNONYMOUS_EXON | python $SCRIPTS/filter_bed.py > "$safe_tmp_dir/exon.bed"

            $BEDTOOLS/bin/bedtools subtract -a "$safe_tmp_dir/intron.bed" -b "$safe_tmp_dir/exon.bed" > $output.bed

            rm -r "$safe_tmp_dir"
        """

        branch.COMBINED_SYNONYMOUS = output.bed
    }
}

build_capture_stats = {
    stage_status("build_capture_stats", "enter", "n/a");
    output.dir = "qc"
    produce( "exon_coverage_stats.txt" ) {
        exec """
            python $SCRIPTS/calculate_exon_coverage.py --capture $EXOME_TARGET --exons $BASE/designs/genelists/exons.bed > qc/exon_coverage_stats.txt
        """
    }
    stage_status("build_capture_stats", "exit", "n/a");
}

generate_pipeline_id = {
    doc "Generate a pipeline run ID for this batch"
    output.dir="results"
    produce("run_id") {
      exec """
        python $SCRIPTS/update_pipeline_run_id.py --id $ID_FILE --increment True > $output
      """
    }
   // This line is necessary on some distributed file systems (e.g. MCRI) to ensure that
   // files get synced between nodes
   file("results").listFiles()
   run_id = new File('results/run_id').text.trim()
}

set_target_info = {

    doc "Validate and set information about the target region to be processed"

    var HG19_CHROM_INFO : false

    exec """
        echo "set_target_info: combined target is $COMBINED_TARGET"
    """

    branch.splice_region_window=false
    branch.splice_region_bed_flag=""
    branch.multi_annovar=false

    branch.batch = batch
    branch.target_name = branch.name // this is the analysis profile
    branch.target_bed_file = "../design/${target_name}.bed"
    branch.target_gene_file = "../design/${target_name}.genes.txt"
    branch.target_samples = sample_info.grep { it.value.target == target_name }*.value*.sample
    branch.transcripts_file = "../design/${target_name}.transcripts.txt"
    branch.target_config = "../design/${target_name}.settings.txt"

    println "Checking for target gene file: $target_gene_file"
    produce(target_gene_file) {
        exec """
            cp $BASE/designs/$target_name/${target_name}.genes.txt $target_gene_file;
        """
    }

    println "Checking for target bed file: $target_bed_file"
    produce(target_bed_file) {
        exec """
            if [ -e $BASE/designs/$target_name/${target_name}.bed ];
            then
                python $SCRIPTS/combine_target_regions.py --bedfiles $BASE/designs/$target_name/${target_name}.bed --genefiles_required ../design/${target_name}.addonce.*.genes.txt --exons $BASE/designs/genelists/exons.bed > $target_bed_file;
            else
                python $SCRIPTS/combine_target_regions.py --genefiles $target_gene_file --genefiles_required ../design/${target_name}.addonce.*.genes.txt --exons $BASE/designs/genelists/exons.bed > $target_bed_file;
            fi
        """
    }

    produce(transcripts_file) {
        exec """
            if [ -e $BASE/designs/$target_name/${target_name}.transcripts.txt ];
            then
                cp $BASE/designs/$target_name/${target_name}.transcripts.txt $transcripts_file;
            else
                touch $transcripts_file;
            fi
        """
    }

    produce(target_config) {
        exec """
            if [ -e $BASE/designs/$target_name/${target_name}.settings.txt ];
            then
                cp $BASE/designs/$target_name/${target_name}.settings.txt $output.txt;
            else
                touch $output.txt;
            fi
        """
    }

    exec """
        if [ -e $BASE/designs/$target_name/${target_name}.pgx.vcf ] && [ ! -e ../designs/${target_name}.pgx.vcf ];
        then
            cp $BASE/designs/$target_name/${target_name}.pgx.vcf ../designs;
        fi
    """

    // if we do this, don't use "using", due to a bug in bpipe https://github.com/ssadedin/bpipe/issues/179
    // Load arbitrary settings related to the target
    println "Loading settings for target region $branch.name from ${file(target_config).absolutePath}"
    load file(target_config).absolutePath

    if(branch.multi_annovar) {
        println "Enabling multiple Annovar annotation sources for $target_name"
        branch.annovar = multiple_annovar 
    }
    
    println "Target $target_name is processing samples $target_samples"
}

create_splice_site_bed = {
    // note: this stage requires annovar

    // If no splice region window is defined, simply set the 
    if(!splice_region_window) {
        branch.splice_region_bed_flag = ""
        return
    }

    msg "Setting regions for calling / annotation of splice variants to $splice_region_window bp past exon boundary"

    output.dir="../design"
    produce(target_name + ".splice.bed", target_name + ".exons.bed") {
        exec """
            python $SCRIPTS/create_exon_bed.py  
                -c -s $target_bed_file $ANNOVAR_DB/hg19_refGene.txt $transcripts_file -
              | $BEDTOOLS/bin/bedtools slop -g $HG19_CHROM_INFO -b $splice_region_window -i - > $output.bed

            python $SCRIPTS/create_exon_bed.py  
                -c $target_bed_file $ANNOVAR_DB/hg19_refGene.txt $transcripts_file $output2.bed
        """

        branch.splice_region_bed_flag = "-L $output1.bed"
        branch.exon_bed_file = output2.bed
    }

    println "Splice region flat = $splice_region_bed_flag"
    println "Exon bed file = $exon_bed_file"
}

sample_similarity_report = {

    doc "Create a report indicating the difference in count of variants for each combination of samples"

    output.dir = "qc"

    produce("similarity_report.txt") {
        exec """
            $JAVA -Xmx4g -cp "${GROOVY_HOME}/embeddable/groovy-all-${GROOVY_VERSION}.jar:$BASE/tools/java_libs/*" VCFSimilarity $inputs.vcf > $output.txt
             """
    }
}

validate_batch = {
    doc "Validates batch results"
    String diseaseGeneLists = ANALYSIS_PROFILES.collect { "$BASE/designs/${it}/${it}.genes.txt" }.join(" ")
    produce("results/missing_from_exons.genes.txt", "results/${run_id}_batch_validation.md", "results/${run_id}_batch_validation.html") {
      exec """
          cat ../design/*.genes.txt | python $SCRIPTS/find_missing_genes.py $BASE/designs/genelists/exons.bed > results/missing_from_exons.genes.txt

          if [ -e $BASE/designs/genelists/annovar.bed ]; then
            cat ../design/*.genes.txt | python $SCRIPTS/find_missing_genes.py $BASE/designs/genelists/annovar.bed > results/missing_from_annovar.genes.txt;
          fi

          if [ -e $BASE/designs/genelists/incidentalome.genes.txt ]; then
            python $SCRIPTS/validate_genelists.py --exclude $BASE/designs/genelists/incidentalome.genes.txt $diseaseGeneLists > results/excluded_genes_analyzed.txt;
          fi

          python $SCRIPTS/validate_batch.py --missing_exons results/missing_from_exons.genes.txt --missing_annovar results/missing_from_annovar.genes.txt --excluded_genes results/excluded_genes_analyzed.txt > results/${run_id}_batch_validation.md

          python $SCRIPTS/markdown2.py --extras tables < results/${run_id}_batch_validation.md | python $SCRIPTS/prettify_markdown.py > results/${run_id}_batch_validation.html
      """, "validate_batch"
    }
}

write_run_info = {
    doc "write out all versions that are relevant to this particular run"
    output.dir = "results"

    produce("${run_id}_pipeline_run_info.txt") {
        exec """
            python $SCRIPTS/write_run_info.py --run_id ${run_id} --base "$BASE" > $output.txt
        """
    }
}

create_sample_metadata = {
    doc "Create a new samples.txt file that includes the pipeline ID"
    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    output.dir="results"
    produce("results/samples.meta") {
        from(sample_metadata_file) {
            exec """
                python $SCRIPTS/update_pipeline_run_id.py --id results/run_id --parse True < $sample_metadata_file > results/samples.meta
            """
        }
    }
}

generate_ped_files = {
    doc "Generate necessary PED files from the details in the sample metadata file"
    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    output.dir="results"

    produce("${run_id}_families.log") {
        from(sample_metadata_file) {
            exec """
                mkdir -p results

                python $SCRIPTS/generate_peds.py --prefix "results/${run_id}_family_" < $sample_metadata_file 2>&1 | tee "results/${run_id}_families.log"
            """
        }
    }
}

update_sample_database = {
    doc "Write details of this sample analysis to the database (if specified)"

    stage_status("update_sample_database", "enter", sample);

    if (SAMPLE_DB && SAMPLE_DB != '') {
        exec """
            python $SCRIPTS/update_sample_db.py --db "$SAMPLE_DB" --sample "${sample}" --run_id "${run_id}" --analysis "${analysis}" --capture "${EXOME_TARGET}" --pipeline_version "`cat $BASE/version.txt`"
        """, "update_sample_database"
    }
    else {
        stage_status("update_sample_database", "skipping...", sample);
    }

    stage_status("update_sample_database", "exit", sample);
}

mark_batch_finished = {
    stage_status("mark_batch_finished", "enter", "");

    var commandline: "";

    if (POST_ANALYSIS_READ_ONLY == true) {
        commandline += "--read_only ";
    }
    if (POST_ANALYSIS_MOVE == true) {
        commandline += "--move ";
    }
    
    exec """
        python $SCRIPTS/mark_batch_finished.py $commandline
    """
    stage_status("mark_batch_finished", "exit", "");
}

///////////////////////////////////////////////////////////////////
// segments
///////////////////////////////////////////////////////////////////

initialize_batch_run = segment {
    // ANALYSIS_PROFILES = sample_info*.value*.target as Set
    // Check the basic sample information first
    check_sample_info +  // check that fastq files are present
    check_tools +
    update_gene_lists + // build new gene lists by adding sample specific genes to cohort

    // Create a single BED that contains all the regions we want to call variants in
    create_combined_target + 
    create_synonymous_target + // regions where synonymous snvs are not filtered
    build_capture_stats + // how well covered genes are by the capture

    generate_pipeline_id + // make a new pipeline run ID file if required

    generate_ped_files // generate ped files required for trio analysis
}

finish_batch_run = segment {
   // report on similarity between samples
   sample_similarity_report +

   // check overall quality of results
   validate_batch +

   // write all genelist versions to results
   write_run_info +

   // update metadata and pipeline ID
   create_sample_metadata +

   // move and mark as read only
   mark_batch_finished
}

// configure target regions and other settings for each profile
initialize_profiles = segment {
    set_target_info // + 
    // create_splice_site_bed // currently disabled, not used (and requires annovar)
}

