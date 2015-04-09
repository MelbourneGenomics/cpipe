// vim: ts=4:expandtab:sw=4:cindent
//
// A pipeline for creating a set of FASTQ reads that contain spiked in
// mutations in a number of clinically interesting ways
//
// NOTE: requires patched version of bamsurgeon that accepts 6 columns of input!!!!!

requires BASE : "Please set the BASE variable to the root of the Cpipe distribution"

load "$BASE/pipeline/config.groovy"
load "$BASE/pipeline/pipeline_stages_config.groovy"

slop_bed = {
    exec """
        $BEDTOOLS/bin/bedtools slop -b 100 -g $REFBASE/hg19.genome -i $BASE/designs/CARDIOM/CARDIOM.bed > $output.bed
     """

    branch.extract_region = output.bed
}

extract_regions = {
    exec """
        $SAMTOOLS/samtools view -b $input.bam -L $input.bed > $output.bam
     """
}

add_variants = {
    requires VARIANT_BED : "BED-ish file specifying variants to spike into sample"
    exec """

        export PATH=`dirname $BWA`:"$PATH"

        cut -f 1,2,3,4,5 $VARIANT_BED | $PYTHON $BAMSURGEON/addsnv.py -v /dev/stdin --bamfile $input.bam -r $REF -o ${output.bam.prefix}.tmp.bam
       
        $SAMTOOLS/samtools sort ${output.bam.prefix}.tmp.bam $output.bam.prefix
        
        rm -rf *_logs_*.bam
    """
}


sam_to_fastq = {

    requires extract_region : "A BED file specifying the region to extract from the BAM file"

    from(extract_region) transform(".bam") to("_R1.fastq.gz", "_R2.fastq.gz") {
        exec """
            $SAMTOOLS/samtools view -b $input.bam -L $input.bed |
            java -Xmx2g -jar $PICARD_HOME/lib/SamToFastq.jar
                    TMP_DIR=$TMP_DIR
                    I=/dev/stdin
                    F=${output1.prefix}
                    F2=${output2.prefix}
                    VALIDATION_STRINGENCY=LENIENT || true

            gzip $output1.prefix

            gzip $output2.prefix

        """
    }
}

run { 
    slop_bed + extract_regions + index_bam + add_variants + index_bam + sam_to_fastq
}
