#!/bin/bash
# parameters
# * 1: input
# * 2: output
# * 3: htslib
# * 4: vep
# * 5: tools
# * 6: condel
# * 7: dbnsfp

VARIANTS=`grep -c -v '^#' < $1`
echo "$VARIANTS variant(s) found in $1"
if [ $VARIANTS -eq 0 ];
then
    grep '^#' $1 > $2
else
    PATH="$PATH:$3"
    perl $4/variant_effect_predictor.pl \
        --allele_number \
        --cache \
        --canonical \
        --check_alleles \
        --check_existing \
        --dir $4/../vep_cache \
        --dir_plugins $5/vep_plugins \
        --fasta $4/../vep_cache/homo_sapiens/83_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        --force_overwrite \
        --gmaf \
        --hgvs \
        -i $1 \
        --maf_1kg \
        --maf_esp \
        --maf_exac \
        -o $2 \
        --offline \
        --per_gene \
        --plugin Condel,$6/config,s \
        --plugin dbNSFP,$7/dbNSFP.gz,LRT_score,LRT_pred,phyloP7way_vertebrate,GERP++_RS,CADD_raw,CADD_phred,MutationTaster_score,MutationTaster_pred,1000Gp3_AC,1000Gp3_AF,ExAC_AC,ExAC_AF,clinvar_rs,clinvar_clnsig,clinvar_trait,Polyphen2_HDIV_score,Polyphen2_HDIV_pred \
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
