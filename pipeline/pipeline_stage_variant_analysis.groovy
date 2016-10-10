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
        python $SCRIPTS/filter_tsv.py --ad ${HARD_FILTER_AD} --af ${HARD_FILTER_AF} --dp ${HARD_FILTER_DP} --qual ${HARD_FILTER_QUAL} --proband ${sample} < $input > $output
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
    doc "Annotate variants using VEP"

    def DBNSFP_OPTS=""
    if(ENABLE_DBNSFP) {
        DBNSFP_OPTS="--plugin dbNSFP,$DBNSFP/dbNSFP.gz,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,VEST3_score,VEST3_rankscore,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,LRT_Omega,UniSNP_ids,1000Gp1_AC,1000Gp1_AF,1000Gp1_AFR_AC,1000Gp1_AFR_AF,1000Gp1_EUR_AC,1000Gp1_EUR_AF,1000Gp1_AMR_AC,1000Gp1_AMR_AF,1000Gp1_ASN_AC,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF ARIC5606_AA_AC,ARIC5606_AA_AF,ARIC5606_EA_AC,ARIC5606_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,clinvar_rs,clinvar_clnsig,clinvar_trait,COSMIC_ID,COSMIC_CNT"
    }

    stage_status("vcf_annotate", "enter", "${sample} ${branch.analysis}");
    output.dir="variants"
    
    // NB we write an empty output file for the case where vep doesn't write anything
    //    /usr/bin/env bash $SCRIPTS/vcf_annotate.sh "$input.vcf" "$output.vcf" "$HTSLIB" "$VEP" "$TOOLS" "$CONDEL" "$DBNSFP"
    exec """
        export PERL5LIB="$PERL5LIB:$TOOLS/perl5:$TOOLS/perl5/lib/perl5";

        echo "$VARIANTS variant(s) found in $input.vcf";
        VARIANTS=`grep -c -v '^#' < $input.vcf`;
        if [ $VARIANTS -eq 0 ] ;
        then
            grep '^#' $input.vcf > $output.vcf;
        else
            PATH="$PATH:$HTSLIB";
            perl $VEP/variant_effect_predictor.pl
                --allele_number
                --assembly GRCh37
                --cache
                --canonical
                --check_alleles
                --check_existing
                --dir $VEP_CACHE
                --dir_plugins $TOOLS/vep_plugins
                --fasta $VEP_CACHE/homo_sapiens_refseq/85_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
                --force_overwrite
                --gmaf
                --hgvs
                -i $input.vcf
                --maf_1kg
                --maf_esp
                --maf_exac
                -o $output.vcf
                --offline
                --plugin Condel,$CONDEL/config,s ${DBNSFP_OPTS}
                --plugin Grantham
                --polyphen b
                --protein
                --pubmed
                --refseq
                --sift b
                -species homo_sapiens
                --symbol
                --vcf
                --vcf_info_field ANN
                --verbose;
        fi
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
    // /usr/bin/env bash $SCRIPTS/vcf_post_annotation_filter.sh "$input.vcf" "$output.vcf" "$VEP" "$TOOLS"
    exec """
        VARIANTS=`grep -c -v '^#' < $input.vcf`

        echo "$VARIANTS variant(s) found in $input.vcf"

        if [ $VARIANTS -eq 0 ];
        then
          cp $input.vcf $output.vcf;
        else
          PERL5LIB="$TOOLS/perl5:$TOOLS/perl5/lib/perl5:$VEP" perl $VEP/filter_vep.pl 
            --input_file $input.vcf
            --filter "Consequence not matches stream" 
            --filter "BIOTYPE match protein_coding"
            --filter "Feature" 
            --force_overwrite
            --format vcf
            -o $output.vcf
            --only_matched ;
        fi
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
