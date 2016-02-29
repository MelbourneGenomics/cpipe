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
                   -dcov 1600 
                   -l INFO 
                   -L $COMBINED_TARGET $splice_region_bed_flag
                   --interval_padding $INTERVAL_PADDING_CALL
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
                   --interval_padding $INTERVAL_PADDING_CALL
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

//////////////////////////////////////////////////////////////////////
// segments
//////////////////////////////////////////////////////////////////////

variant_discovery = segment {
   call_variants_trio + 
   call_pgx + 
   merge_pgx +
   filter_variants + 
   merge_variants
}

