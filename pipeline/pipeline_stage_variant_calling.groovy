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
                   -L $COMBINED_TARGET 
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
   
    // old implementation
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
                   -L $COMBINED_TARGET 
                   --interval_padding $INTERVAL_PADDING_CALL
                   -A AlleleBalance -A Coverage -A FisherStrand 
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}

//////////////////////////////////////////////////////////////////
call_variants_individual_vcf = {
    doc "Call variants and output to VCF"
    stage_status("call_variants_individual_vcf", "enter", sample); 
    output.dir="variants"

    // haplotype calling when not part of a trio
    // java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --emitRefConfidence GVCF --num_cpu_threads_per_data_thread 1 -A AlleleBalance -A GCContent -A GenotypeSummaries -A LikelihoodRankSumTest -A StrandBiasBySample -A VariantType -I $sample\.merge.dedup.realign.recal.bam -L /vlsci/VR0320/shared/production/1.0.4/designs/nextera_rapid_capture_exome_1.2/target_regions.bed -o $sample\.hap.raw.g.vcf --bamOutput $sample\.HC.bam -log $sample\.log -ped txxxx.ped --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf

    var call_conf:5.0, 
        emit_conf:5.0

    // transform("bam") to ("hc.vcf", "hc.bam") {
    produce("${sample}.hc.vcf", "${sample}.hc.bam") {
        from("bam") {
            exec """
                $JAVA -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller 
                    -R $REF 
                    --num_cpu_threads_per_data_thread 1 
                    -A AlleleBalance 
                    -A GCContent 
                    -A GenotypeSummaries 
                    -A LikelihoodRankSumTest 
                    -A StrandBiasBySample 
                    -A VariantType 
                    -I $input.bam
                    -L $COMBINED_TARGET 
                    -o $output.hc.vcf 
                    --bamOutput $output.hc.bam 
                    --logging_level INFO
                    --dbsnp $DBSNP 
                    --interval_padding $INTERVAL_PADDING_CALL
                    -stand_call_conf $call_conf 
                    -stand_emit_conf $emit_conf
            """, "gatk_call_variants"
        }
    }
    stage_status("call_variants_individual_vcf", "exit", sample); 
}

//////////////////////////////////////////////////////////////////
call_variants_individual_gvcf = {
    doc "Call variants and output to GVCF"
    stage_status("call_variants_individual_gvcf", "enter", sample); 
    output.dir="variants"

    // haplotype calling when not part of a trio
    //java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --emitRefConfidence GVCF --num_cpu_threads_per_data_thread 1 -G Standard -A AlleleBalance -A AlleleBalanceBySample -A DepthPerAlleleBySample -A GCContent -A GenotypeSummaries -A HardyWeinberg -A HomopolymerRun -A LikelihoodRankSumTest -A LowMQ -A MappingQualityZero -A SampleList -A SpanningDeletions -A StrandBiasBySample -A TandemRepeatAnnotator -A VariantType -A TransmissionDisequilibriumTest -I $sample\.merge.dedup.realign.recal.bam -L /vlsci/VR0320/shared/production/1.0.4/designs/nextera_rapid_capture_exome_1.2/target_regions.bed -o $sample\.hap.raw.g.vcf --bamOutput $sample\.HC.bam -log $sample\.log -ped txxxx.ped --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf
    // 23-may-2016
    // java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --emitRefConfidence GVCF --num_cpu_threads_per_data_thread 1 -A AlleleBalance -A GCContent -A GenotypeSummaries -A LikelihoodRankSumTest -A StrandBiasBySample -A VariantType -I $sample\.merge.dedup.realign.recal.bam -L /vlsci/VR0320/shared/production/1.0.4/designs/nextera_rapid_capture_exome_1.2/target_regions.bed -o $sample\.hap.raw.g.vcf --bamOutput $sample\.HC.bam -log $sample\.log -ped txxxx.ped --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf

    // transform("bam") to ("hc.g.vcf", "hc.bam") {
    produce("${sample}.hc.g.vcf", "${sample}.hc.bam") {
        from("bam") {
            exec """
                $JAVA -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller 
                    -R $REF 
                    --emitRefConfidence GVCF 
                    --num_cpu_threads_per_data_thread 1 
                    -A AlleleBalance 
                    -A GCContent 
                    -A GenotypeSummaries 
                    -A LikelihoodRankSumTest 
                    -A StrandBiasBySample 
                    -A VariantType 
                    -I $input.bam
                    -L $COMBINED_TARGET 
                    -o $output.hc.g.vcf 
                    --bamOutput $output.hc.bam 
                    --logging_level INFO
                    --dbsnp $DBSNP 
                    --interval_padding $INTERVAL_PADDING_CALL
            """, "gatk_call_variants"
        }
    }
    stage_status("call_variants_individual_gvcf", "exit", sample); 
}

//////////////////////////////////////////////////////////////////
call_variants_trio_gvcf = {
    doc "Call variants with a PED file and output to GVCF"
    output.dir="variants"

    // first stage of trio analysis => includes the ped file
    // TODO is this separation required?

    // java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --emitRefConfidence GVCF --num_cpu_threads_per_data_thread 1 -G Standard -A AlleleBalance -A AlleleBalanceBySample -A DepthPerAlleleBySample -A GCContent -A GenotypeSummaries -A HardyWeinberg -A HomopolymerRun -A LikelihoodRankSumTest -A LowMQ -A MappingQualityZero -A SampleList -A SpanningDeletions -A StrandBiasBySample -A TandemRepeatAnnotator -A VariantType -A TransmissionDisequilibriumTest -I $sample\.merge.dedup.realign.recal.bam -L /vlsci/VR0320/shared/production/1.0.4/designs/nextera_rapid_capture_exome_1.2/target_regions.bed -o $sample\.hap.raw.gvcf --bamOutput $sample\.HC.bam -log $sample\.log -ped txxxx.ped --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf
    // 23-may-2016
    // java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --emitRefConfidence GVCF --num_cpu_threads_per_data_thread 1 -A AlleleBalance -A GCContent -A GenotypeSummaries -A LikelihoodRankSumTest -A StrandBiasBySample -A VariantType -I $sample\.merge.dedup.realign.recal.bam -L /vlsci/VR0320/shared/production/1.0.4/designs/nextera_rapid_capture_exome_1.2/target_regions.bed -o $sample\.hap.raw.g.vcf --bamOutput $sample\.HC.bam -log $sample\.log -ped txxxx.ped --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf
    transform("bam", "ped") to ("hc.g.vcf", "hc.bam") {
        exec """
            $JAVA -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller 
                -R $REF 
                --emitRefConfidence GVCF 
                --num_cpu_threads_per_data_thread 1 
                -A AlleleBalance 
                -A GCContent 
                -A GenotypeSummaries 
                -A LikelihoodRankSumTest 
                -A StrandBiasBySample 
                -A VariantType 
                -I $input.bam
                -L $COMBINED_TARGET 
                -o $output.hc.g.vcf 
                --bamOutput $output.hc.bam 
                --logging_level INFO
                -ped $input.ped 
                --dbsnp $DBSNP
                --interval_padding $INTERVAL_PADDING_CALL
        """, "gatk_call_variants"
    }
}

// For the legacy GATK, we must fall back to UnifiedGenotyper
// for calling variants
call_variants_gatk = GATK_LEGACY ? call_variants_ug : call_variants_hc

call_pgx = {
    doc "Call Pharmacogenomic variants using GATK Unified Genotyper"
    output.dir="variants"

    stage_status("call_pgx", "enter", sample)

    var call_conf:5.0, 
        emit_conf:5.0

    if(!file("../design/${target_name}.pgx.vcf").exists()) {
        stage_status("call_pgx", "exit (returning)", sample)
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
    stage_status("call_pgx", "exit", sample)
}

//////////////////////////////////////////////////////////////////////
merge_pgx = {
    doc "Merge multiple VCF files into one file"
    output.dir="variants"

    stage_status("merge_pgx", "enter", target_name)

    if(!file("../design/${target_name}.pgx.vcf").exists()) {
        stage_status("merge_pgx", "forwarding", target_name)
        if (GATK_VARIANT_ONLY) {
            forward input.hc.vcf.toString() // workaround for Bpipe bug
        }
        else {
            forward input.hc.g.vcf.toString() // workaround for Bpipe bug
        }
        stage_status("merge_pgx", "exit (returning)", target_name)
        return
    }

    msg "Merging vcf files: " + inputs.vcf
    if (GATK_VARIANT_ONLY) {
        exec """
            $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T CombineVariants
            -R $REF
            --variant $input.hc.vcf
            --variant $input.pgx.vcf
            --out $output.vcf
         """
    }
    else {
        exec """
            $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T CombineVariants
            -R $REF
            --variant $input.hc.g.vcf
            --variant $input.pgx.vcf
            --out $output.vcf
         """
    }

    stage_status("merge_pgx", "exit", target_name)
}

//////////////////////////////////////////////////////////////////////
// segments
//////////////////////////////////////////////////////////////////////

if (GATK_VARIANT_ONLY) {
    variant_discovery = segment {
        call_variants_individual_vcf + 
        call_pgx + 
        merge_pgx
    }
}
else {
    variant_discovery = segment {
        call_variants_individual_gvcf + 
        call_pgx + 
        merge_pgx
    }
}
