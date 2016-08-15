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
    
    // generate a custom bed file that only includes the incidentalome for this sample
    def sample_bed_file = "$target_bed_file.${sample}.bed"
    produce(sample_bed_file) {
        exec """
            python $SCRIPTS/combine_target_regions.py --genefiles $target_gene_file --genefiles_required ../design/${target_name}.addonce.${sample}.genes.txt --exons $BASE/designs/genelists/exons.bed --bedfiles $BASE/designs/${target_name}/${target_name}.bed > $output.bed
        """
    }
 
    println "Processing input files ${files} for target region $sample_bed_file"

    forward files
}

fastqc = {
    doc "Run FASTQC to generate QC metrics for raw reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip') {
        exec "$FASTQC/fastqc --extract -o ${output.dir} $inputs.gz"
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
           cat fastqc/${sample}_*fastqc/summary.txt |
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
    stage_status("align_bwa", "enter", sample);

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

    // var BWA_THREADS: false;

    if(!BWA_THREADS) {
        BWA_THREADS = 1
    }

    stage_status("align_bwa", "output file is ${outputFile}", sample);
    produce(outputFile) {
        
        uses(threads:BWA_THREADS) {
            //    Note: the results are filtered with flag 0x100 because bwa mem includes multiple 
            //    secondary alignments for each read, which upsets downstream tools such as 
            //    GATK and Picard.
            def safe_tmp_dir = [TMPDIR, UUID.randomUUID().toString()].join( File.separator )
            exec """
                    set -o pipefail
    
                    mkdir "$safe_tmp_dir"
    
                    $BWA mem -M -t $threads -k $seed_length 
                             -R "@RG\\tID:${sample}_${lane}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample_info[sample].library}\\tSM:${sample}"  
                             $REF $input1.gz $input2.gz | 
                             $SAMTOOLS/samtools view -F 0x100 -bSu - | $SAMTOOLS/samtools sort -o ${output.prefix}.bam -T "$safe_tmp_dir/bamsort"
    
                    rm -r "$safe_tmp_dir"
            ""","bwamem"
        }
    }
    stage_status("align_bwa", "exit", sample);
}

index_bam = {
    stage_status("index_bam", "enter", sample);
 
    doc "Create an index for a BAM file"
 
    // A bit of a hack to ensure the index appears in the
    // same directory as the input bam, no matter where it is
    // nb: fixed in new version of Bpipe
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
    transform("bam") to ("bam.bai") {
        exec "$SAMTOOLS/samtools index $input.bam"
    }
    stage_status("index_bam", "forwarding", sample);
    forward input
    stage_status("index_bam", "exit", sample);
}
 
merge_bams = {
    doc """
        Merge the BAM files from multiple lanes together.
        <p>
        When there is only one BAM file the merge is skipped, and a 
        simple copy is made instead
        """

    stage_status("merge_bams", "enter", sample);
    output.dir="align"

    produce(sample + ".merge.bam") {
        // If there is only 1 bam file, then there is no need to merge,
        // just alias the name 
        //if(inputs.bam.size()==1)  {
        //    alias(input.bam) to(output.bam)
            // msg "Skipping merge of $inputs.bam because there is only one file"
            // This use of symbolic links may be questionable
            // However if the ordinary case involves only one
            // bam file then there may be some significant savings
            // from doing this.
            // exec "ln -sf ${file(input.bam).name} $output.bam; ln -sf ${file(input.bam).name}.bai ${output.bam}.bai;"
        //}
        //else {
            msg "Merging $inputs.bam size=${inputs.bam.size()}"
            exec """
                $JAVA -Xmx2g -jar $PICARD_HOME/picard.jar MergeSamFiles
                    ${inputs.bam.withFlag("INPUT=")}
                    VALIDATION_STRINGENCY=LENIENT
                    ASSUME_SORTED=true
                    CREATE_INDEX=true
                    OUTPUT=$output.bam
             """, "merge"
        //}
    }
    stage_status("merge_bams", "exit", sample);
}

dedup = {
    doc "Remove PCR duplicates from reads"
    stage_status("dedup", "enter", sample);
    output.dir="align"

    var MAX_DUPLICATION_RATE : 30

    def safe_tmp_dir = [TMPDIR, UUID.randomUUID().toString()].join( File.separator )
    exec """
        mkdir -p "$safe_tmp_dir"

        $JAVA -Xmx4g -Djava.io.tmpdir=$safe_tmp_dir -jar $PICARD_HOME/picard.jar MarkDuplicates
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
    stage_status("dedup", "exit", sample);
}

cleanup_initial_bams = {
    cleanup("*.merge.bam", ~".*L00[0-9].bam")
}

cleanup_intermediate_bams = {
    cleanup("*.dedup.bam", "*.realign.bam")
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
    stage_status("recal", "enter", sample);
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
    stage_status("recal", "exit", sample);
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
            $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
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
                $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
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
        recal_count + 
        recal
    }
}

///////////////////////////////////////////////////////////////////
// segments
///////////////////////////////////////////////////////////////////
// all samples do this
analysis_ready_reads = segment {
    set_sample_info +
    "%.gz" * [ fastqc ] + check_fastqc +
     ~"(.*)_R[0-9][_.].*fastq.gz" * 
     [ 
         trim_fastq + 
         align_bwa + index_bam + 
         cleanup_trim_fastq 
     ] +
     merge_bams +
     dedup + 
     // cleanup_initial_bams + // seems to mess things up
     realignIntervals + 
     realign + index_bam +
     bsqr_recalibration + index_bam
     // cleanup_intermediate_bams
}

