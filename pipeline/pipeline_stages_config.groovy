// vim: ts=4:sw=4:expandtab:cindent:number
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

ENABLE_CADD=true

call_variants_ug = {
    doc "Call SNPs/SNVs using GATK Unified Genotyper"
    output.dir="variants"

    // Default values of quality thresholds
    var call_conf:5.0, 
        emit_conf:5.0

    transform("bam","bam") to("metrics.txt","vcf") {
        exec """
                $JAVA -Xmx6g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
                   -R $REF 
                   -I $input.bam 
                   -nt 4
                   --dbsnp $DBSNP 
                   -stand_call_conf $call_conf -stand_emit_conf $emit_conf
                   -dcov 1600 
                   -l INFO 
                   -L $COMBINED_TARGET $splice_region_bed_flag
                   -A AlleleBalance -A FisherStrand 
                   -glm BOTH
                   -metrics $output.txt
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}



set_target_info = {

    doc "Validate and set information about the target region to be processed"

    var HG19_CHROM_INFO : false

    branch.splice_region_window=false
    branch.splice_region_bed_flag=""
    branch.multi_annovar=false

    branch.batch = batch
    branch.target_name = branch.name
    branch.target_bed_file = "../design/${target_name}.bed"
    branch.target_gene_file = "../design/${target_name}.genes.txt"
    branch.target_samples = sample_info.grep { it.value.target == target_name }*.value*.sample
    branch.transcripts_file = "../design/${target_name}.transcripts.txt"
    branch.target_config = "../design/${target_name}.settings.txt"

    println "Checking for target bed file : $target_bed_file"

    produce(target_bed_file) {
        exec """
                cp $BASE/designs/$target_name/${target_name}.bed $target_bed_file; 
        """
    }

    produce(target_gene_file) {
        exec """
            cp $BASE/designs/$target_name/${target_name}.genes.txt $target_gene_file;
        """
    }

    produce(transcripts_file) {
        exec """
            cp $BASE/designs/$target_name/${target_name}.transcripts.txt $transcripts_file;
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

    // Load arbitrary settings related to the target
    println "Loading settings for target region $branch.name from ${file(target_config).absolutePath}"
    load file(target_config).absolutePath
    if(multi_annovar) {
        println "Enabling multiple Annovar annotation sources for $target_name"
        branch.annoar = multiple_annovar 
    }
    
    println "Target $target_name is processing samples $target_samples"
}

init_analysis_profile = {
  // This stage is a placeholder to allow individual analysis profiles
  // to perform initialization steps
}

create_splice_site_bed = {

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

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    if(sample_info[sample].target != target_name) {
        // This is expected because every file is processed for every target/flagship
        succeed "No files to process for sample $sample, target $target_name"
    }

    // Patient specific variants are not supported yet
    // If they are provided, we should not process the patient at all
    check {
        if(sample_info[sample].variantsFile?.trim()) 
                exec "false" // force failure
    } otherwise { 
        succeed """
             Study $sample is configured with a sample specific variant file. The pipeline currently does 
             not support sample specific variants. Please remove the variant file from the configuration
             to allow processing.
         """.trim().stripIndent() to channel: cpipe_operator, subject: "Invalid configuration for Study $sample"
    }

    def files = sample_info[sample].files.fastq

    println "Processing input files ${files} for target region ${target_bed_file}"
    forward files
}

check_tools = {
    doc """
        Checks for presence of optional tools and sets appropriate pipeline variables
        to enable or disable corresponding pipeline features
        """

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
}

check_sample_info = {

    doc "Validate basic sample information is correct"

    def missingSummary = []
    for(sample in samples) {

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
}

create_combined_target = {

    // Construct the region for variant calling from 
    //
    //   a) all of the disease cohort BED files
    //   b) the EXOME target regions
    //
    // This way we avoid calling variants over the entire genome, but still
    // include everything of interest
    String diseaseBeds = ANALYSIS_PROFILES.collect{"$BASE/designs/${it}/${it}.bed"}.join(",")

    output.dir = "../design"

    produce("combined_target_regions.bed") {
        exec """
            cat $diseaseBeds $EXOME_TARGET | 
                cut -f 1,2,3 | 
                $BEDTOOLS/bin/bedtools sort | 
                $BEDTOOLS/bin/bedtools merge > $output.bed
        """

        branch.COMBINED_TARGET = output.bed
    }
}

fastqc = {
    doc "Run FASTQC to generate QC metrics for raw reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip') {
        exec "$FASTQC/fastqc -o ${output.dir} $inputs.gz"
    }
}

check_fastqc = {

    doc "Search for any failures in FastQC output and abort further processing if they are found"

    check {
       // NOTE: we remove per-base-sequence content and
       // per-base-gc-content from examination because Nextera
       // appears to contain natural biases that flag QC failures 
       // here.
       exec """
           cat fastqc/"${sample}"_*fastqc/summary.txt |
               grep -v "Per base sequence content" |
               grep -v "Per base GC content" |
               grep -q 'FAIL' && exit 1

           exit 0
       ""","local"
    } otherwise {
        if(CHECK_FASTQC_FAILURES) {
            succeed report('templates/fastqc_failure.html') to channel: cpipe_operator, 
                                                            subject: "Sample $sample has failed FastQC Check", 
                                                            file: input.zip
        } else {
            send report('templates/fastqc_failure.html') to channel: cpipe_operator, 
                                                            subject: "Sample $sample has failed FastQC Check", 
                                                            file: input.zip
        }
    }

    check("FASTQ Format") {
        exec """
            awk -F'\\t' '/Illumina/ { where=match(\$2, /[0-9.]+/); { result=substr(\$2, where, RLENGTH); exit(result<1.7); } }' fastqc/${sample}_*_fastqc/fastqc_data.txt
        ""","local"
    } otherwise {
        println "=" * 100
        println "Sample $sample is encoded using a quality encoding incompatible with this pipeline."
        println "Please convert the data first using maq ill2sanger."
        println "=" * 100

        succeed report('templates/fastqc_failure.html') to channel: cpipe_operator, 
                                                        subject: "Sample $sample is encoded with incompatible quality scores (Illumina < 1.7)", 
                                                        file: input.zip
    }
}

trim_fastq = {
   output.dir="align"
   if(ADAPTERS_FASTA) {
       filter("trim","trim") {
            exec """
                $JAVA -Xmx4g -jar $TRIMMOMATIC/trimmomatic-0.30.jar PE -phred33 
                    $input1.gz $input2.gz 
                    $output1.gz ${output1.prefix}.unpaired.gz 
                    $output2.gz ${output2.prefix}.unpaired.gz 
                    ILLUMINACLIP:$ADAPTERS_FASTA:2:40:15 LEADING:3 TRAILING:6 SLIDINGWINDOW:4:15 MINLEN:36
            """
       }
   }
}

cleanup_trim_fastq = {
    if(ADAPTERS_FASTA)
        cleanup "*.fastq.trim.gz"
}

align_bwa = {

    doc "Align with bwa mem algorithm."

    output.dir = "align"

    var seed_length : 19

    def lanes = inputs.gz.collect { (it.toString() =~ /_(L[0-9]{1,3})_/)[0][1] }.unique()
    if(lanes.size()!=1) 
        succeed report('templates/invalid_input.html') to channel: cpipe_operator, 
                                                          subject: "Invalid input files for sample $sample: Bad lane information",
                                                          message: """Failed to identify a unique lane number from FASTQ files: ${inputs.gz}. 
                                                                        Please check the format of the input file names""".stripIndent()
        
    branch.lane = lanes[0]

    def outputFile = sample + "_" + Hash.sha1(inputs.gz*.toString().join(",")) + "_" + lane + ".bam"
    produce(outputFile) {
        //    Note: the results are filtered with flag 0x100 because bwa mem includes multiple 
        //    secondary alignments for each read, which upsets downstream tools such as 
        //    GATK and Picard.
        exec """
                set -o pipefail

                $BWA mem -M -t $threads -k $seed_length 
                         -R "@RG\\tID:${sample}_${lane}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample_info[sample].library}\\tSM:${sample}"  
                         $REF $input1.gz $input2.gz | 
                         $SAMTOOLS/samtools view -F 0x100 -bSu - | $SAMTOOLS/samtools sort - $output.prefix
        ""","bwamem"
    }
}

merge_pgx = {
    doc "Merge multiple VCF files into one file"
    output.dir="variants"

    if(!file("../design/${target_name}.pgx.vcf").exists()) {
        forward input.recal.vcf.toString() // workaround for Bpipe bug
        return
    }

    msg "Merging vcf files: " + inputs.vcf
    exec """
            $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T CombineVariants
            -R $REF
            --variant $input.recal.vcf
            --variant $input.pgx.vcf
            --out $output.vcf
         """
}

merge_vcf = {
    doc "Merge multiple VCF files into one file"
    output.dir="variants"

    produce(target_name + ".merge.vcf") {
        msg "Merging vcf files: " + inputs.vcf
        exec """
                $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
                -T CombineVariants
                -R $REF
                ${inputs.vcf.withFlag("--variant")}
                --out $output.vcf
             """
    }
}

merge_target_vcfs = {

    doc "Merges only the VCFs belonging to the current analysis profile"

    output.dir="variants"

    var enable_family_excel : false
    if(!enable_family_excel)
        succeed "Family VCF output not configured for $target_name"

    // Find the sample names that are in the current analysis profile
    def target_samples = sample_info.grep { it.value.target == target_name }*.value*.sample

    // Work with raw filtered VCFs, not ones that were already annotated
    def variant_files = target_samples*.plus(".*.filter.vcf")

    from(variant_files) produce(target_name + ".merge.vcf") {
        msg "Merging vcf files: " + inputs.vcf
        exec """
                $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
                -T CombineVariants
                -R $REF
                ${inputs.vcf.withFlag("--variant")}
                --out $output.vcf
             """
    }
}

merge_bams = {
    doc """
        Merge the BAM files from multiple lanes together.
        <p>
        When there is only one BAM file the merge is skipped, and a 
        simple copy is made instead
        """

    output.dir="align"

    produce(sample + ".merge.bam") {
        // If there is only 1 bam file, then there is no need to merge,
        // just alias the name 
        if(inputs.bam.size()==1)  {
           // alias(input.bam) to(output.bam)
            msg "Skipping merge of $inputs.bam because there is only one file"
            // This use of symbolic links may be questionable
            // However if the ordinary case involves only one
            // bam file then there may be some significant savings
            // from doing this.
            exec "ln -s ${file(input.bam).name} $output.bam; ln -s ${file(input.bam).name}.bai ${output.bam}.bai;"
        }
        else {
            msg "Merging $inputs.bam size=${inputs.bam.size()}"
            exec """
                $JAVA -Xmx2g -jar $PICARD_HOME/lib/MergeSamFiles.jar
                    ${inputs.bam.withFlag("INPUT=")}
                    VALIDATION_STRINGENCY=LENIENT
                    ASSUME_SORTED=true
                    CREATE_INDEX=true
                    OUTPUT=$output.bam
             """, "merge"
        }
    }
}

index_bam = {

    doc "Create an index for a BAM file"

    // A bit of a hack to ensure the index appears in the
    // same directory as the input bam, no matter where it is
    // nb: fixed in new version of Bpipe
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
    transform("bam") to ("bam.bai") {
        exec "$SAMTOOLS/samtools index $input.bam"
    }
    forward input
}

realignIntervals = {
    doc "Discover candidate regions for realignment in an alignment with GATK"
    output.dir="align"
    exec """
        $JAVA -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
            -T RealignerTargetCreator 
            -R $REF 
            -I $input.bam 
            -L $COMBINED_TARGET $splice_region_bed_flag
            --known $GOLD_STANDARD_INDELS 
            -o $output.intervals
    """, "realign_target_creator"
}

realign = {
    doc "Apply GATK local realignment to specified intervals in an alignment"
    output.dir="align"
    exec """
        $JAVA -Xmx5g -jar $GATK/GenomeAnalysisTK.jar 
             -T IndelRealigner 
             -R $REF 
             -I $input.bam 
             -L $COMBINED_TARGET $splice_region_bed_flag
             -targetIntervals $input.intervals 
             -o $output.bam
    ""","local_realign"
}

dedup = {
    doc "Remove PCR duplicates from reads"
    output.dir="align"

    var MAX_DUPLICATION_RATE : 30

    def safe_tmp_dir = [TMPDIR, UUID.randomUUID().toString()].join( File.separator )
    exec """
        mkdir -p "$safe_tmp_dir"

        $JAVA -Xmx4g -Djava.io.tmpdir=$safe_tmp_dir -jar $PICARD_HOME/lib/MarkDuplicates.jar
             INPUT=$input.bam 
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$output.metrics
             CREATE_INDEX=true
             OUTPUT=$output.bam

        rm -r "$safe_tmp_dir"
    """

    check {
        exec """
            DUPLICATION_RATE=`grep -A 1 LIBRARY $output.metrics | cut -f 8 | tail -1  | awk '{ print int(\$1 * 100) }'`

            [ $DUPLICATION_RATE -lt $MAX_DUPLICATION_RATE ]

        ""","local"
    } otherwise {
        send text {"Rate of PCR duplicates for sample $sample is higher than $MAX_DUPLICATION_RATE"} to channel: cpipe_operator 
    }
}

recal_count = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="align"
    INDEL_QUALS=""
    // To use lite version of GATK uncomment below
    // INDEL_QUALS="--disable_indel_quals"
    exec """
        $JAVA -Xmx5g -jar $GATK/GenomeAnalysisTK.jar 
             -T BaseRecalibrator 
             -I $input.bam 
             -R $REF 
             -L $COMBINED_TARGET $splice_region_bed_flag
             --knownSites $DBSNP $INDEL_QUALS
             -l INFO 
             -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate 
             -o $output.counts
    """, "recalibrate_bam"
}

recal = {
    doc "Apply recalibration quality adjustments so that quality scores match actual observed error rates"
    output.dir="align"
    exec """
          $JAVA -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
               -T PrintReads 
               -I $input.bam 
               -BQSR $input.counts 
               -L $COMBINED_TARGET $splice_region_bed_flag
               -R $REF 
               -l INFO 
               -o $output.bam
        """, "recalibrate_bam"
}

legacy_recal_count = {
    doc title: "Calculate factors correlated with poor base quality for use in recalibrating reads",
        desc: """Includes recalibration using the following standard covariates:
                    ContextCovariate
                    CycleCovariate
                    QualityScoreCovariate
                    ReadGroupCovariate
             """,
        inputs: "BAM file containing reads, VCF file containing known variants from dbSNP (132)",
        outputs: "CSV file containing correlation information"

    output.dir="align"

    transform("recal.csv") {
        exec """
            java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T BaseRecalibrator
            -R $REF
            -l INFO
            -L $COMBINED_TARGET $splice_region_bed_flag
            -I $input.bam
            --disable_indel_quals
            -knownSites $DBSNP
            -o $output
            ""","count_covariates"
    }
}

legacy_recal = {
    msg "Performing recalibration ..."
    output.dir="align"
    from("csv","bam") {
        transform('bam') {
            exec """
                java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
                    -l INFO
                    -L $COMBINED_TARGET $splice_region_bed_flag
                    -R $REF
                    -I $input.bam
                    -T PrintReads
                    -BQSR $input.csv
                    -o $output.bam
                ""","recalibrate_bam"
        }
    }
}

if(GATK_LEGACY) {
    bsqr_recalibration = segment {
        legacy_recal_count + legacy_recal
    }
}
else {
    bsqr_recalibration = segment {
        recal_count + recal
    }
}

cleanup_initial_bams = {
    cleanup("*.merge.bam", ~".*L00[0-9].bam")
}

cleanup_intermediate_bams = {
    cleanup("*.dedup.bam", "*.realign.bam")
}

call_variants_hc = {
    doc "Call SNPs/SNVs using GATK Unified Genotyper"
    output.dir="variants"

    // Default values of confidence thresholds
    // come from the Broad web site. However
    // these may be higher than suitable in our context
    var call_conf:5.0, 
        emit_conf:5.0

    transform("bam") to("vcf") {
        exec """

            $JAVA -Xmx6g -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller
                   -R $REF 
                   -I $input.bam 
                   --dbsnp $DBSNP 
                   -stand_call_conf $call_conf -stand_emit_conf $emit_conf
                   -dcov 1600 
                   -l INFO 
                   -L $COMBINED_TARGET $splice_region_bed_flag
                   -A AlleleBalance -A Coverage -A FisherStrand 
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}

// For the legacy GATK, we must fall back to UnifiedGenotyper
// for calling variants
call_variants_gatk = GATK_LEGACY ? call_variants_ug : call_variants_hc

call_pgx = {
    doc "Call Pharmacogenomic variants using GATK Unified Genotyper"
    output.dir="variants"

    var call_conf:5.0, 
        emit_conf:5.0

    if(!file("../design/${target_name}.pgx.vcf").exists())
        return

    transform("bam","bam") to("metrics","pgx.vcf") {
        exec """
                $JAVA -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
                   -R $REF 
                   -I $input.bam 
                   -nt 4
                   --dbsnp $DBSNP 
                   -stand_call_conf $call_conf -stand_emit_conf $emit_conf
                   --output_mode EMIT_ALL_SITES
                   -dcov 1600 
                   -l INFO 
                   -L ../design/${target_name}.pgx.vcf
                   -A AlleleBalance -A Coverage -A FisherStrand 
                   -glm BOTH
                   -metrics $output.metrics
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}

filter_variants = {
    doc "Select only variants in the genomic regions defined for the $target_name target"
    output.dir="variants"

    def pgx_flag = ""
    if(file("../design/${target_name}.pgx.vcf").exists()) {
        pgx_flag = "-L ../design/${target_name}.pgx.vcf"
    }

    filter("filter") {
        exec """
            $JAVA -Xmx2g -jar $GATK/GenomeAnalysisTK.jar 
                 -R $REF
                 -T SelectVariants 
                 --variant $input.vcf 
                 -L $target_bed_file $splice_region_bed_flag $pgx_flag
                 -o $output.vcf 
        """
    }
}

@filter("vep")
annotate_vep = {
    doc "Annotate variants using VEP to add Ensemble annotations"
    output.dir="variants"

    // Note: if the input VCF file is empty, VEP will not create an output.
    // To avoid this causing the pipeline to fail, we first copy only the
    // headers from the input file to the output file so that if VEP does not
    // overwrite the file, we end up with an empty file
    exec """
        grep '^#' $input.vcf > $output.vcf 

        PERL5LIB="$CONDEL:\$PERL5LIB"
        perl $VEP/variant_effect_predictor.pl --cache --dir $VEP/../vep_cache 
            -i $input.vcf 
            --vcf -o $output.vcf 
            -species human 
            --canonical --per_gene --protein 
            --sift=b --polyphen=b
            --symbol hgnc --force_overwrite --hgvs  --maf_1kg --maf_esp --pubmed
            --plugin Condel,$CONDEL/config,s
    """, "vep"
}

@filter("snpeff")
annotate_snpeff = {
    output.dir="variants"

    var enable_snpeff:false,
        SNPEFF : false

    if(!enable_snpeff)
        succeed "Snpeff support not enabled"

    if(!SNPEFF)
        fail "Please define the SNPEFF variable to point to the location of your SNPEFF installation"

    exec """
            $JAVA -Xmx2g -jar $SNPEFF/snpEff.jar eff -c $SNPEFF/snpEff.config -treatAllAsProteinCoding false -a 2 hg19 $input.vcf  > $output.vcf
    ""","snpeff"
}

calc_coverage_stats = {
    doc "Calculate coverage across a target region using Bedtools"
    output.dir="qc"

    var MIN_ONTARGET_PERCENTAGE : 50

    transform("bam","bam") to(file(target_bed_file).name+".cov.gz","ontarget.txt") {
        exec """
          $BEDTOOLS/bin/coverageBed -d  -abam $input.bam -b $target_bed_file | gzip -c > $output.gz

          $SAMTOOLS/samtools view -L $COMBINED_TARGET $input.bam | wc | awk '{ print \$1 }' > $output2.txt
        """
    }
}

check_ontarget_perc = {
    var MIN_ONTARGET_PERCENTAGE : 50
    check {
        exec """
            RAW_READ_COUNT=`cat $input.ontarget.txt`

            ONTARGET_PERC=`grep -A 1 LIBRARY $input.metrics | tail -1 | awk '{ print int(((\$3 * 2) / $RAW_READ_COUNT))*100 }'`

            [ $ONTARGET_PERC -lt $MIN_ONTARGET_PERCENTAGE ]

             """
    } otherwise {
        send text {"On target read percentage for $sample < $MIN_ONTARGET_PERCENTAGE"} to channel: cpipe_operator 
    }
}

check_coverage = {

    output.dir = "qc"

    def medianCov
    transform("cov.gz") to("cov.stats.median", "cov.stats.csv") {

        R {"""
            bam.cov = read.table(pipe("gunzip -c $input.cov.gz"), col.names=c("chr","start","end", "gene", "offset", "cov"))
            meds = aggregate(bam.cov$cov, list(bam.cov$gene), median)
            write.csv(data.frame(Gene=meds[,1],MedianCov=meds$x), "$output.csv", quote=F, row.names=F)
            writeLines(as.character(median(bam.cov$cov)), "$output.median")
        """}

        // HACK to ensure file sync on distributed file system
        file(output.dir).listFiles()
        medianCov = Math.round(file(output.median).text.toFloat())
    }

    check {
        exec "[ $medianCov -ge $MEDIAN_COVERAGE_THRESHOLD ]"
    } otherwise {
        // It may seem odd to call this a success, but what we mean by it is that
        // Bpipe should not fail the whole pipeline, merely this branch of it
        succeed report('templates/sample_failure.html') to channel: cpipe_operator, 
                                                           median: medianCov, 
                                                           file:output.csv, 
                                                           subject:"Sample $sample has failed with insufficient median coverage ($medianCov)"
    }
}

check_karyotype = {

    doc "Compare the inferred sex of the sample to the inferred karyotype from the sequencing data"

    def karyotype_file = "results/" + run_id + '_' + sample + '.summary.karyotype.tsv'
    check {
        exec """
            [ `grep '^Sex' $karyotype_file | cut -f 2` == "UNKNOWN" ] || [ `grep '^Sex' $karyotype_file | cut -f 2` == `grep 'Inferred Sex' $karyotype_file | cut -f 2` ]
        """
    } otherwise {
        // It may seem odd to call this a success, but what we mean by it is that
        // Bpipe should not fail the whole pipeline, merely this branch of it
        succeed report('templates/sample_failure.html') to channel: cpipe_operator, 
                                                           median: medianCov, 
                                                           file: karyotype_file,
                                                           subject:"Sample $sample has a different sex than inferred from sequencing data"
     }
}

augment_condel = {

    doc "Extract Condel scores from VEP annotated VCF files and add them to Annovar formatted CSV output"

    output.dir="variants"
    from("*.hg19_multianno*.csv") filter("con") {
        exec """
            JAVA_OPTS="-Xmx4g -Djava.awt.headless=true" $GROOVY 
                -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar 
                $SCRIPTS/merge_condel.groovy
                    -i $input.vcf
                    -a $input.csv
                    -o $output.csv
        """
    }
}

index_vcf = {
    output.dir="variants"
    var sort_vcf : true
    if(sort_vcf) {
        transform("vcf") to("sort.vcf") {
            exec """

                $IGVTOOLS/igvtools sort $input.vcf $output.vcf 

                if [ ! -e $output.vcf ];
                then
                    grep '^#' $input.vcf > $output.vcf ;
                fi

                $IGVTOOLS/igvtools index $output.vcf
            """
        }
    }
    else {
        transform("vcf") to("vcf.idx") {
            exec """
                $IGVTOOLS/igvtools index $input.vcf
            """
        }
    }
}

vcf_to_excel = {

    doc "Convert a VCF output file to Excel format, merging information from Annovar"

    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)",
             ANNOTATION_VARIANT_DB : "File name of SQLite variant database for storing variants"

    var exclude_variant_types : "synonymous SNV",
        out_of_cohort_filter_threshold : OUT_OF_COHORT_VARIANT_COUNT_FILTER

    check {
        exec "ls results/${target_name}.qc.xlsx > /dev/null 2>&1"
    } otherwise { 
        succeed "No samples succeeded for target $target_name" 
    }

    def pgx_flag = ""
    if(file("../design/${target_name}.pgx.vcf").exists()) {
        pgx_flag = "-pgx ../design/${target_name}.pgx.vcf"
    }

    // Default behavior: exclude all synonymous variants from output
    def EXCLUDE_VARIANT_TYPES=exclude_variant_types.toString().split(",")*.trim().join(",")

    println "Excluding variant types: $EXCLUDE_VARIANT_TYPES"
    println "Filtering out variants observed more than $out_of_cohort_filter_threshold times"

    output.dir="results"

    def all_outputs = [target_name + ".xlsx"] + target_samples.collect { run_id + '_' + it + ".annovarx.csv" }
    from("*.hg19_multianno.*.csv", "*.vcf") produce(all_outputs) {
        exec """
            echo "Creating $outputs.csv"

            JAVA_OPTS="-Xmx12g -Djava.awt.headless=true" $GROOVY 
                -cp $SCRIPTS:$GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/vcf_to_excel.annovar.groovy 
                -s '${target_samples.join(",")}'
                ${inputs.csv.withFlag("-a")}
                ${inputs.vcf.withFlag("-vcf")}
                -x "$EXCLUDE_VARIANT_TYPES"
                -db $ANNOTATION_VARIANT_DB
                -o $output.xlsx
                -oocf $out_of_cohort_filter_threshold
                -si $sample_metadata_file
                -gc $target_gene_file ${pgx_flag}
                -annox $output.dir
                -log ${target_name}_filtering.log
                -prefix $run_id
                ${inputs.bam.withFlag("-bam")}
        """, "vcf_to_excel"
    }
}

vcf_to_family_excel = {

    doc """
         Convert variants annotated with SnpEff and Annovar to an excel based format
         designed for diagnostics from family based sequencing.
        """

    requires ANNOTATION_VARIANT_DB : "File name of SQLite variant database for storing variants"

    output.dir = "results"

    def UNIQUE = unique.toBoolean() ? " -unique " : ""

    def target_samples = sample_info.grep { it.value.target == target_name }*.value*.sample
    
    def annovar_files = target_samples*.plus(".*.hg19_multianno.con.sig.csv") 

    println "Annovar files for target $target_name are " + annovar_files
    from(annovar_files) {
        produce(target_name + ".family.xlsx") {
            exec """
                JAVA_OPTS="-Xmx6g -Djava.awt.headless=true" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/vcf_to_excel.family.groovy 
                    -p "" 
                    -db $ANNOTATION_VARIANT_DB
                    -ped $input.ped
                    -o $output.xlsx
                    -p $run_id
                    $UNIQUE $input.vcf $inputs.csv 
            """, "vcf_to_family_excel"
        }
    }
}

family_vcf = segment { merge_target_vcfs + annotate_snpeff + index_vcf.using(sort_vcf:false) + vcf_to_family_excel }

plot_coverage = {
    doc "Create plots showing coverage distributions for alignment"
    from("cov.gz") {
        transform("cum.png") {
            msg "Plotting coverage for sample $sample_name"
            R {"""
                sample = '$input.gz'
                name = '$sample'
                samplecov = read.table(file=pipe(paste("gunzip -c ",sample," | awk '{ print $NF }' ", sep='')))
                tf = table(factor(samplecov$V1, levels=0:max(samplecov)))
                cs = cumsum(tf)
                png("$output.png")
                plot(1-cs/max(cs), ylim=c(0,1), type="l", xaxs="i", xlim=c(0,800),
                     xaxt="n",
                     main=paste("Cum. Coverage Distribution for ", name), xlab='Coverage Depth', ylab='Frequency')
                axis(side=1, at=seq(0,800,25))
                dev.off()
            """} 
        }
    }
}

gatk_depth_of_coverage = {

    doc "Calculate statistics about depth of coverage for an alignment using GATK"

    output.dir = "qc"
    transform("bam") to("."+target_name + ".cov.sample_cumulative_coverage_proportions", 
                         "."+target_name + ".cov.sample_interval_statistics") { 
        exec """
            $JAVA -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
               -R $REF
               -T DepthOfCoverage 
               -o $output.sample_cumulative_coverage_proportions.prefix
               --omitDepthOutputAtEachBase
               -I $input.bam
               -ct 1 -ct 10 -ct 20 -ct 50 -ct 100
               -L $target_bed_file
        """
    }
}

insert_size_metrics = {

    doc "Generates statistics about distribution of DNA fragment sizes"

    var MIN_MEDIAN_INSERT_SIZE : 70,
        MAX_MEDIAN_INSERT_SIZE : 240

    output.dir="qc"
    exec """
        $JAVA -Xmx4g -jar $PICARD_HOME/lib/CollectInsertSizeMetrics.jar INPUT=$input.bam O=$output.txt H=$output.pdf
    """

    check {
        exec """
              INSERT_SIZE=`grep -A 1 MEDIAN_INSERT_SIZE $output.txt | cut -f 1 | tail -1`

              echo "Median insert size = $INSERT_SIZE"

              [ $INSERT_SIZE -gt $MIN_MEDIAN_INSERT_SIZE ] && [ $INSERT_SIZE -lt $MAX_MEDIAN_INSERT_SIZE ]
             """, "local"
    } otherwise {
        send text {"""
            WARNING: Insert size distribution for $sample has median out of 
            range $MIN_MEDIAN_INSERT_SIZE - $MAX_MEDIAN_INSERT_SIZE
        """} to channel: cpipe_operator, file: output.pdf
    }
}

qc_excel_report = {

    doc "Create an excel file containing a summary of QC data for all the samples for a given target region"

    var LOW_COVERAGE_THRESHOLD : 15,
        LOW_COVERAGE_WIDTH : 1

    output.dir="results"

    def samples = sample_info.grep { it.value.target == target_name }.collect { it.value.sample }
    produce(target_name + ".qc.xlsx") {
            exec """
                JAVA_OPTS="-Xmx16g -Djava.awt.headless=true" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/qc_excel_report.groovy 
                    -s ${target_samples.join(",")} 
                    -t $LOW_COVERAGE_THRESHOLD
                    -w $LOW_COVERAGE_WIDTH
                    -low qc ${inputs.dedup.metrics.withFlag('-metrics')}
                    -o $output.xlsx
                    $inputs.sample_cumulative_coverage_proportions  
                    $inputs.sample_interval_statistics 
                    $inputs.gz
            ""","qc_excel_report"
    }
}

@filter("sig")
annotate_significance = {
    doc "Add clinical significance category annotations as defined by Melbourne Genomics"
        var MAF_THRESHOLD_RARE : 0.01,
            MAF_THRESHOLD_VERY_RARE : 0.0005,
            CONDEL_THRESHOLD : 0.7

        output.dir="variants"
        from("con.csv") {
            exec """
                python $SCRIPTS/annotate_significance.py 
                -a $input.csv
                -f $MAF_THRESHOLD_RARE
                -r $MAF_THRESHOLD_VERY_RARE
                -c $CONDEL_THRESHOLD
                > $output.csv
                """
        }
}

annovar_table = {

    output.dir="variants"

    transform("vcf","vcf") to("av", "hg19_multianno.csv") {
        exec """
            $ANNOVAR/convert2annovar.pl $input.vcf -format vcf4 > $output.av

            $ANNOVAR/table_annovar.pl $output.av $ANNOVAR_DB/  -buildver hg19 
            -protocol refGene,phastConsElements46way,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,exac03,snp138,ljb26_all
            -operation g,r,r,f,f,f,f,f 
            -nastring . 
            --otherinfo   
            --csvout
            --outfile $output.csv.prefix.prefix
            --argument '-exonicsplicing -splicing $splice_region_window',,,,,,,

            sed -i '/^Chr,/ s/\\.refGene//g' $output.csv
        """
    }
}

merge_annovar_reports = {
    msg "Augmenting Annovar output with extra columns ($input.csv) ..."
    msg "Inputs are $inputs"
    output.dir="variants"
    from(["*.refgene.*.csv", "*.knowngene.*.csv"]) filter("merge") {
        exec "python $SCRIPTS/merge_knowngene_annotations.py $input1 $input2 > $output.csv"
    }
}



add_to_database = {

    doc "Add discovered variants to a database to enable annotation of number of observations of the variant"

    requires UPDATE_VARIANT_DB : "SQLite database to store variants in"

    output.dir="variants"

    uses(variantdb:1) {
        exec """

            echo "====> Adding variants for flaship $target_name to database"

            JAVA_OPTS="-Xmx24g" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/vcf_to_db.groovy 
                   -v $input.recal.vcf 
                   -a $input.csv 
                   -db $UPDATE_VARIANT_DB 
                   -cohort $target_name
                   -idmask '$SAMPLE_ID_MASK'
                   -b "$batch"

            echo "<==== Finished adding variants for flaship $target_name to database"

            echo "Variants from $input.recal.vcf were added to database $VARIANT_DB on ${new Date()}" > $output.txt
        """, "add_to_database"
    }
}

reorder = {
    filter('reorder') {
        exec """
            $JAVA -Xmx2g -jar $PICARD_HOME/lib/ReorderSam.jar
                TMP_DIR=$TMP_DIR
                I=$input.bam
                O=$output.bam
                VALIDATION_STRINGENCY=LENIENT
                REFERENCE=$REF
            """ 
    }
}

summary_pdf = {

    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    output.dir="results"

    produce("${run_id}_${sample}.summary.pdf","${run_id}_${sample}.summary.karyotype.tsv") {

        // -metrics $input.metrics
        exec """
             JAVA_OPTS="-Xmx3g" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar $SCRIPTS/qc_pdf.groovy 
                -cov $input.cov.gz
                -ontarget $input.ontarget.txt ${inputs.metrics.withFlag("-metrics")}
                -study $sample 
                -meta $sample_metadata_file
                -threshold 20 
                -classes GOOD:95:GREEN,PASS:80:ORANGE,FAIL:0:RED 
                -exome $EXOME_TARGET
                -gc $target_gene_file ${pgx_flag}
                -bam $input.bam
                -o $output.pdf
        """

        branch.karyotype = output.tsv

        send text {"Sequencing Results for Study $sample"} to channel: cpipe_operator, file: output.pdf
    }
}

exon_qc_report = {

    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    output.dir="results"

    var enable_exon_report : false

    if(!enable_exon_report)  {
        msg "Exon level coverage report not enabled for $target_name"
        return
    }

    produce("${sample}.exon.qc.xlsx", "${sample}.exon.qc.tsv") {
        exec """
             JAVA_OPTS="-Xmx3g" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/exon_qc_report.groovy 
                -cov $input.cov.gz
                -targets $target_bed_file
                -refgene $ANNOVAR_DB/hg19_refGene.txt 
                -x $output.xlsx
                -o $output.tsv
        """
    }
}

sample_similarity_report = {

    doc "Create a report indicating the difference in count of variants for each combination of samples"

    output.dir = "qc"

    produce("similarity_report.txt") {
        exec """
            $JAVA -Xmx4g -cp $GROOVY_HOME/embeddable/groovy-all-2.3.4.jar:$GROOVY_NGS/groovy-ngs-utils.jar VCFSimilarity $inputs.vcf > $output.txt
             """
    }
}

provenance_report = {
    branch.sample = branch.name
    output.dir = "results"
    produce(run_id + '_' + sample + ".provenance.pdf") {
       send report("scripts/provenance_report.groovy") to file: output.pdf
    }
}

annovar_to_lovd = {
    branch.sample = branch.name
    output.dir="results/lovd"
    produce(run_id + '_' + sample +"_LOVD") {
        exec """
            python $SCRIPTS/annovar2LOVD.py --csv $input.annovarx.csv --meta $sample_metadata_file --dir results/lovd
        """
    }
}

generate_pipeline_id = {
    doc "Generate a pipeline run ID for this batch"
    output.dir="results"
    produce("results/run_id") {
      exec """
        python $SCRIPTS/update_pipeline_run_id.py --id $ID_FILE --increment True > results/run_id
      """
      run_id = new File('results/run_id').text.trim()
    }
}

create_sample_metadata = {
    doc "Create a new samples.txt file that includes the pipeline ID"
    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    output.dir="results"
    produce("results/samples.meta") {
      exec """
          python $SCRIPTS/update_pipeline_run_id.py --id results/run_id --parse True < $sample_metadata_file > results/samples.meta
      """
    }
}

variant_bams = {

    doc "Create a bam file for each variant containing only reads overlapping 100bp either side of that variant"

    output.dir = "results/variant_bams"

    from(branch.name+'*annovarx.csv', branch.name+'.*.recal.bam') {   
        // Slight hack here. Produce a log file that bpipe can track to confirm that the bams were produced.
        // Bpipe is not actually tracking the variant bams themselves. 
        produce(branch.name + ".variant_bams_log.txt") {
            exec """
                python $SCRIPTS/variant_bams.py --bam $input.bam --csv $input.csv --outdir $output.dir --log $output.txt --samtoolsdir $SAMTOOLS
            """
        }
    }
}

@filter("aug")
augment_transcript_ids = {
        doc "Add an additional column indicating the VCGS transcript identifier for the VCGS transcript / isoform affected by each variant"
                output.dir="variants"

                    msg "Augmenting Annovar output with extra columns ($input.csv) ..."
                        exec "python $SCRIPTS/augment_transcripts.py $transcripts_file $input.csv $input.exonic_variant_function > $output.csv"
}

/*
   multiple_annovar = segment {
       [
               annovar_summarize.using(source:"refgene") + add_splice_variants + augment_transcript_ids,
                       annovar_summarize.using(source:"knowngene") + add_splice_variants + augment_transcript_ids
       ] + merge_annovar_reports
}
*/


