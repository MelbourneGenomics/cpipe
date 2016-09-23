#!/bin/bash
# parameters
# * 1: input
# * 2: output
# * 3: htslib
# * 4: vep
# * 5: tools
# * 6: condel
# * 7: dbnsfp

# 23-may-2016 perl /vlsci/VR0320/shared/production/2.2.0/tools/vep/83/variant_effect_predictor.pl 
# --allele_number 
# --cache 
# --canonical 
# --check_existing 
# --check_alleles  
# --dir /vlsci/VR0320/shared/production/2.2.0/tools/vep/vep_cache 
# --dir_plugins /vlsci/VR0320/shared/production/2.2.0/tools/vep_plugins 
# --fasta /vlsci/VR0320/shared/production/2.2.0/tools/vep/vep_cache/homo_sapiens/83_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa 
# --force_overwrite 
# --gmaf 
# --hgvs 
# -i $trio\.genotype.raw.split.norm.ChildOnly.vcf 
# --maf_1kg 
# --maf_esp 
# --maf_exac 
# -o $trio\.genotype.raw.split.norm.ChildOnly.vep.83.vcf 
# --offline 
# --per_gene # removed 23-sep-2016
# --plugin Condel,/vlsci/VR0320/shared/production/2.2.0/tools/vep_plugins/condel/2.4/config,s 
# --plugin dbNSFP,/vlsci/VR0320/shared/production/2.2.0/tools/vep_plugins/dbNSFP_2.9.1/dbNSFP.EDIT.gz,LRT_score,LRT_pred,phyloP100way_vertebrate,GERP++_RS,CADD_raw,CADD_phred,MutationTaster_score,MutationTaster_pred,1000Gp1_AC,1000Gp1_AF,ExAC_AC,ExAC_AF,clinvar_rs,clinvar_clnsig,clinvar_trait,Polyphen2_HDIV_score,Polyphen2_HDIV_pred
# --polyphen b 
# --protein 
# --pubmed 
# --refseq 
# --sift b 
# -species homo_sapiens 
# --symbol 
# --vcf 
# --vcf_info_field ANN 

TOOLS="$5"
VARIANTS=`grep -c -v '^#' < $1`

export PERL5LIB="$PERL5LIB:$TOOLS/perl5:$TOOLS/perl5/lib/perl5"

echo "$VARIANTS variant(s) found in $1"
if [ $VARIANTS -eq 0 ];
then
    grep '^#' $1 > $2
else
    PATH="$PATH:$3"
    perl $4/variant_effect_predictor.pl \
        --allele_number \
        --assembly GRCh37 \
        --cache \
        --canonical \
        --check_alleles \
        --check_existing \
        --dir $4/../vep_cache \
        --dir_plugins $TOOLS/vep_plugins \
        --fasta $4/../vep_cache/homo_sapiens/83_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
        --force_overwrite \
        --gmaf \
        --hgvs \
        -i $1 \
        --maf_1kg \
        --maf_esp \
        --maf_exac \
        -o $2 \
        --offline \
        --plugin Condel,$6/config,s \
        --plugin dbNSFP,$7/dbNSFP.gz,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,VEST3_score,VEST3_rankscore,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,LRT_Omega,UniSNP_ids,1000Gp1_AC,1000Gp1_AF,1000Gp1_AFR_AC,1000Gp1_AFR_AF,1000Gp1_EUR_AC,1000Gp1_EUR_AF,1000Gp1_AMR_AC,1000Gp1_AMR_AF,1000Gp1_ASN_AC,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF ARIC5606_AA_AC,ARIC5606_AA_AF,ARIC5606_EA_AC,ARIC5606_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,clinvar_rs,clinvar_clnsig,clinvar_trait,COSMIC_ID,COSMIC_CNT \
        --plugin Grantham \
        --polyphen b \
        --protein \
        --pubmed \
        --refseq \
        --sift b \
        -species homo_sapiens \
        --symbol \
        --vcf \
        --vcf_info_field ANN \
        --verbose
fi
