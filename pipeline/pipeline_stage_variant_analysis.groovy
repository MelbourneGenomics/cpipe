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

@filter("hf")
filter_table = {
    doc "filter tsv based on specified cutoffs"
    stage_status("table_filter", "enter", "${sample} ${branch.analysis}");
    output.dir="results"
    exec """
        python $SCRIPTS/filter_tsv.py --ad ${HARD_FILTER_AD} --af ${HARD_FILTER_AF} --dp ${HARD_FILTER_DP} --qual ${HARD_FILTER_QUAL} < $input > $output
    """
    stage_status("table_filter", "exit", "${sample} ${branch.analysis}");
}

vcf_normalize = {
    doc "split VCF lines so that each line contains one and only one variant, and left-normalize all VCF lines"

    stage_status("vcf_normalize", "enter", "${sample} ${branch.analysis}");
    output.dir="variants"
    produce ("${sample}.${analysis}.genotype.norm.vcf") {
        from("variants/${sample}.${analysis}.refined.vcf") {
            exec """
                $BCFTOOLS/bcftools norm -m -both $input.refined.vcf | $BCFTOOLS/bcftools norm -f $REF - -o $output
            """
        }
    }
    stage_status("vcf_normalize", "exit", "${sample} ${branch.analysis}");
}

vcf_filter_child = {
    doc "remove any hom ref variants"
    stage_status("vcf_filter_child", "enter", "${sample} ${branch.analysis}");
    output.dir="variants"
    // 23-may-2016 java -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T SelectVariants -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --variant txxxx.genotype.raw.split.norm.vcf -select 'vc.getGenotype("00NA12879").isHomVar() || vc.getGenotype("00NA12879").isHet() ' -o txxxx.genotype.raw.split.norm.ChildOnly.vcf
    transform ("norm.vcf") to("soi.vcf") {
        exec """
            $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R $REF --variant $input.norm.vcf -select 'vc.getGenotype("$sample").isHomVar() || vc.getGenotype("$sample").isHet()' -o $output.soi.vcf
        """
    }
    stage_status("vcf_filter_child", "exit", "${sample} ${branch.analysis}");
}

@filter("vep")
vcf_annotate = {
    doc "annotate variants"
    stage_status("vcf_annotate", "enter", "${sample} ${branch.analysis}");
    output.dir="variants"
    // NB we write an empty output file for the case where vep doesn't write anything
    exec """
        /usr/bin/env bash $SCRIPTS/vcf_annotate.sh "$input.vcf" "$output.vcf" "$HTSLIB" "$VEP" "$TOOLS" "$CONDEL" "$DBNSFP"
    """
    stage_status("vcf_annotate", "exit", "${sample} ${branch.analysis}");
}

@filter("post_filter")
vcf_post_annotation_filter = {
    doc "filter the vcf post annotation"
    stage_status("vcf_post_annotation_filter", "enter", "${sample} ${branch.analysis}");
    output.dir="variants"
    // PERL5LIB="/vlsci/VR0320/shared/production/2.2.0/tools/vep/83" perl /vlsci/VR0320/shared/production/2.2.0/tools/vep/83/filter_vep.pl --input_file txxxx.genotype.raw.split.norm.ChildOnly.vep.83.vcf --format vcf --filter "Consequence not matches stream" --only_matched --filter "BIOTYPE match protein_coding" --filter "Feature" -o txxxx.genotype.raw.split.norm.ChildOnly.vep.83.FILTER.vcf
    // first copy input to output, otherwise filter_vep creates an empty file
    exec """
        /usr/bin/env bash $SCRIPTS/vcf_post_annotation_filter.sh "$input.vcf" "$output.vcf" "$VEP"
    """
    stage_status("vcf_post_annotation_filter", "exit", "${sample} ${branch.analysis}");
}

vcf_to_table = {
    doc "convert to tab delimited format"
    stage_status("vcf_to_table", "enter", "${sample} ${branch.analysis}");
    output.dir="variants"
    // java -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T VariantsToTable -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F ABHet -F ABHom -F AC -F AF -F AN -F BaseQRankSum -F DB -F DP -F DS -F Dels -F END -F ExcessHet -F FS -F GC -F GQ_MEAN -F GQ_STDDEV -F HRun -F HW -F HaplotypeScore -F InbreedingCoeff -F LikelihoodRankSum -F LowMQ -F MLEAC -F MLEAF -F MQ -F MQ0 -F MQRankSum -F NCC -F OND -F QD -F RAW_MQ -F RPA -F RU -F ReadPosRankSum -F SOR -F STR -F Samples -F TDT -F VariantType -F ANN -GF AB -GF AD -GF DP -GF GQ -GF GT -GF MIN_DP -GF PGT -GF PID -GF PL -GF RGQ -GF SB -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --allowMissingData --showFiltered -V txxxx.genotype.raw.split.norm.ChildOnly.vep.83.FILTER.vcf -o txxxx.genotype.raw.split.norm.ChildOnly.vep.83.FILTER.table
    exec """
        /usr/bin/env bash $SCRIPTS/vcf_to_table.sh "$input.vcf" "$output.table" "$SCRIPTS" "$JAVA" "$GATK" "$REF"
    """
    stage_status("vcf_to_table", "exit", "${sample} ${branch.analysis}");
}

@filter("acr")
annotate_custom_regions = {
    doc "if a custom bed file is specified, add an annotation marking whether it is in a region"
    stage_status("annotate_custom_regions", "enter", "${sample} ${branch.analysis}");
    output.dir="results"
    if (ANNOTATE_CUSTOM_REGIONS != "") {
        exec """
            python $SCRIPTS/annotate_custom_regions.py --bed $ANNOTATE_CUSTOM_REGIONS < $input > $output
        """
    }
    else { // nothing to do
        exec """
            cp "$input" "$output"
        """
    }
    stage_status("annotate_custom_regions", "exit", "${sample} ${branch.analysis}");
}

table_to_lovd = {
    doc "Explodes out all annotation fields for LOVD compatibility"
    stage_status("table_to_lovd", "enter", "${sample} ${branch.analysis}");
    output.dir="results"
    produce("${run_id}_${sample}.${analysis}.lovd.tsv") {
        exec """
            python $SCRIPTS/convert_to_lovd.py --vcf $input.vcf < $input.table > $output.tsv
        """
    }
    stage_status("table_to_lovd", "exit", "${sample} ${branch.analysis}");
}

variant_analysis = segment {
    vcf_normalize +
    vcf_filter_child +
    vcf_annotate +
    vcf_post_annotation_filter + // vep filter
    vcf_to_table +
    filter_table +
    annotate_custom_regions +
    table_to_lovd
}
