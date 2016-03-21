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

//////////////////////////////////////////////////////////////////////
// stages
//////////////////////////////////////////////////////////////////////
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
                   --interval_padding $INTERVAL_PADDING_CALL
                   -A AlleleBalance -A FisherStrand 
                   -glm BOTH
                   -metrics $output.txt
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}

call_variants_hc = {
    doc "Call SNPs/SNVs using GATK Haplotype Caller"
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
                   -l INFO 
                   -L $COMBINED_TARGET $splice_region_bed_flag
                   --interval_padding $INTERVAL_PADDING_CALL
                   -A AlleleBalance -A Coverage -A FisherStrand 
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}

call_variants_individual_gvcf = {
    doc "Call variants and output to GVCF"
    output.dir="variants"

    // haplotype calling when not part of a trio
    //java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --emitRefConfidence GVCF --num_cpu_threads_per_data_thread 1 -G Standard -A AlleleBalance -A AlleleBalanceBySample -A DepthPerAlleleBySample -A GCContent -A GenotypeSummaries -A HardyWeinberg -A HomopolymerRun -A LikelihoodRankSumTest -A LowMQ -A MappingQualityZero -A SampleList -A SpanningDeletions -A StrandBiasBySample -A TandemRepeatAnnotator -A VariantType -A TransmissionDisequilibriumTest -I $sample\.merge.dedup.realign.recal.bam -L /vlsci/VR0320/shared/production/1.0.4/designs/nextera_rapid_capture_exome_1.2/target_regions.bed -o $sample\.hap.raw.g.vcf --bamOutput $sample\.HC.bam -log $sample\.log -ped txxxx.ped --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf

    transform("bam") to ("g.vcf", "hc.bam") {
        exec """
            java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller 
                -R $REF 
                --emitRefConfidence GVCF 
                --num_cpu_threads_per_data_thread 1 
                -G Standard 
                -A AlleleBalance 
                -A AlleleBalanceBySample 
                -A DepthPerAlleleBySample 
                -A GCContent 
                -A GenotypeSummaries 
                -A HardyWeinberg 
                -A HomopolymerRun 
                -A LikelihoodRankSumTest 
                -A LowMQ 
                -A MappingQualityZero 
                -A SampleList 
                -A SpanningDeletions 
                -A StrandBiasBySample 
                -A TandemRepeatAnnotator 
                -A VariantType 
                -A TransmissionDisequilibriumTest 
                -I $input.bam
                -L $COMBINED_TARGET 
                -o $output.g.vcf 
                --bamOutput $output.hc.bam 
                --logging_level INFO
                --dbsnp $DBSNP
        """, "gatk_call_variants"
    }
}

call_variants_trio_gvcf = {
    doc "Call variants with a PED file and output to GVCF"
    output.dir="variants"

    // first stage of trio analysis

    // java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --emitRefConfidence GVCF --num_cpu_threads_per_data_thread 1 -G Standard -A AlleleBalance -A AlleleBalanceBySample -A DepthPerAlleleBySample -A GCContent -A GenotypeSummaries -A HardyWeinberg -A HomopolymerRun -A LikelihoodRankSumTest -A LowMQ -A MappingQualityZero -A SampleList -A SpanningDeletions -A StrandBiasBySample -A TandemRepeatAnnotator -A VariantType -A TransmissionDisequilibriumTest -I $sample\.merge.dedup.realign.recal.bam -L /vlsci/VR0320/shared/production/1.0.4/designs/nextera_rapid_capture_exome_1.2/target_regions.bed -o $sample\.hap.raw.gvcf --bamOutput $sample\.HC.bam -log $sample\.log -ped txxxx.ped --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf
    transform("bam", "ped") to ("g.vcf", "hc.bam") {
        exec """
            java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller 
                -R $REF 
                --emitRefConfidence GVCF 
                --num_cpu_threads_per_data_thread 1 
                -G Standard 
                -A AlleleBalance 
                -A AlleleBalanceBySample 
                -A DepthPerAlleleBySample 
                -A GCContent 
                -A GenotypeSummaries 
                -A HardyWeinberg 
                -A HomopolymerRun 
                -A LikelihoodRankSumTest 
                -A LowMQ 
                -A MappingQualityZero 
                -A SampleList 
                -A SpanningDeletions 
                -A StrandBiasBySample 
                -A TandemRepeatAnnotator 
                -A VariantType 
                -A TransmissionDisequilibriumTest 
                -I $input.bam
                -L $COMBINED_TARGET 
                -o $output.g.vcf 
                --bamOutput $output.hc.bam 
                --logging_level INFO
                -ped $input.ped 
                --dbsnp $DBSNP
        """, "gatk_call_variants"
    }
}

genotype_likelihoods_individual = {
    doc "Convert an individual GVCF to a vcf"
    // java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --disable_auto_index_creation_and_locking_when_reading_rods --num_threads 1 --variant 00NA12877.hap.raw.g.vcf --variant 00NA12878.hap.raw.g.vcf --variant 00NA12879.hap.raw.g.vcf --out txxxx.genotype.raw.vcf -ped txxxx.ped -log txxxx.GenotypeGVCFs.log --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf -G Standard -A AlleleBalance -A AlleleBalanceBySample -A DepthPerAlleleBySample -A GCContent -A GenotypeSummaries -A HardyWeinberg -A LikelihoodRankSumTest -A MappingQualityZero -A SampleList -A SpanningDeletions -A StrandBiasBySample -A TandemRepeatAnnotator -A VariantType -A TransmissionDisequilibriumTest
    transform("vcf") to ("vcf") {
        exec """
          java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs 
            -R $REF 
            --disable_auto_index_creation_and_locking_when_reading_rods 
            --num_threads 1 
            --variant $input.vcf
            --out $output.vcf
            --logging_level INFO 
            --dbsnp $DBSNP
            -G Standard 
            -A AlleleBalance 
            -A AlleleBalanceBySample 
            -A DepthPerAlleleBySample 
            -A GCContent 
            -A GenotypeSummaries 
            -A HardyWeinberg 
            -A LikelihoodRankSumTest 
            -A MappingQualityZero 
            -A SampleList 
            -A SpanningDeletions 
            -A StrandBiasBySample 
            -A TandemRepeatAnnotator 
            -A VariantType 
            -A TransmissionDisequilibriumTest
        """
    }
}

genotype_likelihoods_trio = {
    doc "Convert the GVCF to a vcf using trio inputs"
    // java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --disable_auto_index_creation_and_locking_when_reading_rods --num_threads 1 --variant 00NA12877.hap.raw.g.vcf --variant 00NA12878.hap.raw.g.vcf --variant 00NA12879.hap.raw.g.vcf --out txxxx.genotype.raw.vcf -ped txxxx.ped -log txxxx.GenotypeGVCFs.log --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf -G Standard -A AlleleBalance -A AlleleBalanceBySample -A DepthPerAlleleBySample -A GCContent -A GenotypeSummaries -A HardyWeinberg -A LikelihoodRankSumTest -A MappingQualityZero -A SampleList -A SpanningDeletions -A StrandBiasBySample -A TandemRepeatAnnotator -A VariantType -A TransmissionDisequilibriumTest
    transform("vcf") to ("vcf") {
      exec """
        java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs 
            -R $REF 
            --disable_auto_index_creation_and_locking_when_reading_rods 
            --num_threads 1 
            --variant 00NA12877.hap.raw.g.vcf 
            --variant 00NA12878.hap.raw.g.vcf 
            --variant 00NA12879.hap.raw.g.vcf 
            --out $output.vcf
            -ped $input.ped 
            --logging_level INFO 
            --dbsnp $DBSNP
            -G Standard 
            -A AlleleBalance 
            -A AlleleBalanceBySample 
            -A DepthPerAlleleBySample 
            -A GCContent 
            -A GenotypeSummaries 
            -A HardyWeinberg 
            -A LikelihoodRankSumTest 
            -A MappingQualityZero 
            -A SampleList 
            -A SpanningDeletions 
            -A StrandBiasBySample 
            -A TandemRepeatAnnotator 
            -A VariantType 
            -A TransmissionDisequilibriumTest
      """
    }
}

// For the legacy GATK, we must fall back to UnifiedGenotyper
// for calling variants
call_variants_gatk = GATK_LEGACY ? call_variants_ug : call_variants_hc

call_pgx = {
    doc "Call Pharmacogenomic variants using GATK Unified Genotyper"
    output.dir="variants"

    exec """
        echo "call_pgx: enter"
    """

    var call_conf:5.0, 
        emit_conf:5.0

    if(!file("../design/${target_name}.pgx.vcf").exists()) {
        exec """
            echo "call_pgx: returning"
        """
        return
    }

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
                   --interval_padding $INTERVAL_PADDING_CALL
                   -A AlleleBalance -A Coverage -A FisherStrand 
                   -glm BOTH
                   -metrics $output.metrics
                   -o $output.vcf
            ""","gatk_call_variants"
    }
    exec """
        echo "call_pgx: exit"
    """
}

filter_variants = {
    doc "Select only variants in the genomic regions defined for the $target_name target"
    output.dir="variants"

    exec """
        echo "filter_variants: enter"
    """

    def pgx_flag = ""
    if(file("../design/${target_name}.pgx.vcf").exists()) {
        pgx_flag = "-L ../design/${target_name}.pgx.vcf"
    }

    msg "Filtering variants - finding INDELs"
    exec """
        java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar 
             -R $REF
             -T SelectVariants 
             --variant $input.vcf 
             -L $target_bed_file.${sample}.bed $pgx_flag
             --interval_padding $INTERVAL_PADDING_SNV
             --selectTypeToInclude SNP --selectTypeToInclude MIXED --selectTypeToInclude MNP --selectTypeToInclude SYMBOLIC --selectTypeToInclude NO_VARIATION
             -o $output.snv
    """

    msg "Filtering variants - finding SNVs"
    exec """
        java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar 
             -R $REF
             -T SelectVariants 
             --variant $input.vcf 
             -L $target_bed_file.${sample}.bed $pgx_flag
             --interval_padding $INTERVAL_PADDING_INDEL
             --selectTypeToInclude INDEL
             -o $output.indel
    """
    exec """
        echo "filter_variants: exit"
    """
}

merge_variants_gvcf = {
    doc "Merge SNVs and INDELs"
    output.dir="variants"

    msg "Merging SNVs and INDELs"
    exec """
            java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T CombineGVCFs
            -R $REF
            --variant:indel $input.indel
            --variant:snv $input.snv
            --out $output.vcf
         """
}


merge_variants = {
    doc "Merge SNVs and INDELs"
    output.dir="variants"

    msg "Merging SNVs and INDELs"
    exec """
            java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T CombineVariants
            -R $REF
            --variant:indel $input.indel
            --variant:snv $input.snv
            --out $output.vcf
            --setKey set
            --genotypemergeoption UNSORTED
         """
}

merge_pgx = {
    doc "Merge multiple VCF files into one file"
    output.dir="variants"

    exec """
        echo "merge_pgx: enter ${target_name}"
    """

    if(!file("../design/${target_name}.pgx.vcf").exists()) {
        exec """
            echo "merge_pgx: forwarding..."
        """
        forward input.recal.g.vcf.toString() // workaround for Bpipe bug
        exec """
            echo "merge_pgx: done"
        """
        return
    }

    msg "Merging vcf files: " + inputs.vcf
    exec """
            $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T CombineVariants
            -R $REF
            --variant $input.recal.g.vcf
            --variant $input.pgx.vcf
            --out $output.vcf
         """

    exec """
        echo "merge_pgx: exit"
    """
}

//////////////////////////////////////////////////////////////////////
// segments
//////////////////////////////////////////////////////////////////////

variant_discovery = segment {
    call_variants_individual_gvcf + 
    call_pgx + 
    merge_pgx +
    filter_variants +
    merge_variants_gvcf
}

