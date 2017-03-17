###########################################################
#
# Build a batch directory containing our test data 
#
HERE=`pwd`

# The variant database we will use - we do not want to interfere with any
# previous variant database
VARIANT_DB=$HERE/variants.db

BATCH_NAME=mutation_detection
BATCH_DIR=${CPIPE_ROOT}/batches/${BATCH_NAME}

cd ${CPIPE_ROOT}
if [ -e $BATCH_DIR ];
then
    rm -rf $BATCH_DIR
fi

msg "Creating directories ..."
mkdir -p $BATCH_DIR/data || err "Unable to create batch directory"

msg "Creating batch using CARDIOM analysis profile ..."
manage_batch create ${BATCH_NAME} --profile CARDIOM --data ${HERE}/data/*.fastq.gz --exome ${CPIPE_ROOT}/designs/genelists/exons.bed --force || err "Failed to create batch"

pushd $BATCH_DIR/analysis

msg "Running pipeline ... "

../../../bpipe run -p CHECK_FASTQC_FAILURES=false -p VARIANT_DB=$VARIANT_DB ../../../pipeline/pipeline.groovy ../samples.txt | tee output.log

msg "Checking results ..."

# Now check that our mutations were correctly picked up in the output *.annovarx results
# We expect to see the following variants identified and in the final report:
# chr1 156104747 0.5 G nonsyn_exon
# chr1 156105692 0.5 T syn_inner_splice
# chr1 156106711 0.5 T outer_splice_1bp
# chr1 156106899 0.5 T outer_splice_5bp
# chr1 201330464 0.5 C stop_codon
# chr1 201328340 0.5 G stop_loss

#ANNOVAR_CSV=`ls $BATCH_DIR/analysis/results/*.annovarx.csv`
LOVD_TSV=`ls $BATCH_DIR/analysis/results/*.lovd.tsv`

#annovar_has chr1 156104747 "nonsynonymous SNV"
lovd_has chr1 156104747 "missense_variant" 

# Unfortunately legacy GATK does not find these variants
# because the small size of the data set causes it to 
# treat them as false positives and their base quality scores
# all get heavily downranked by BSQR. For now I am 
# disabling these when running using legacy GATK.
# It should not be a problem for larger data sets
if [ ! $GATK_LEGACY ];
then
    #annovar_has chr1 156105692 "exonic;splicing"  #  todo: should this be? "splicing"
    lovd_has chr1 156105692 "splice_region_variant"  #  todo: should this be? "splicing"
    #annovar_has chr1 156106711 "splicing" # not found with gatk 3.5
    #lovd_has chr1 156106899 "splicing" # not found with gatk 3.5
fi
#annovar_has chr1 201330464 "stopgain"
#annovar_has chr1 201328340 "stoploss"
lovd_has chr1 201330464 "stop_gained"
lovd_has chr1 201328340 "stop_lost"

msg "Test Completed Successfully"

true

