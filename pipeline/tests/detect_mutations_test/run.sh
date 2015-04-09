###########################################################
#
# Build a batch directory containing our test data 
#
HERE=`pwd`

# The variant database we will use - we do not want to interfere with any 
# previous variant database
VARIANT_DB=$HERE/variants.db

BATCH_DIR=$BASE/batches/mutation_detection

cd $BASE
if [ -e $BATCH_DIR ];
then
    rm -rf $BATCH_DIR
fi

msg "Creating directories ..."
mkdir -p $BATCH_DIR/data || err "Unable to create batch directory"

msg "Copying test data ..."
cp -v $HERE/data/*.fastq.gz $BATCH_DIR/data || err "Unable to copy data to batch directory"

msg "Creating batch using CARDIOM analysis profile ..."
./pipeline/scripts/create_batch.sh `basename $BATCH_DIR` CARDIOM || err "Failed to create batch"

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

ANNOVAR_CSV=`ls $BATCH_DIR/analysis/results/*.annovarx.csv`

annovar_has chr1 156104747 "nonsynonymous SNV" 

# This fails because it is classified as synonymous and filtered out
# following up with clinicians to see what the right handling is
annovar_has chr1 156105692 "exonic;splicing"  #  todo: should this be? "splicing"
annovar_has chr1 156106711 "splicing"
annovar_has chr1 156106899 "splicing"
annovar_has chr1 201330464 "stopgain SNV"
annovar_has chr1 201328340 "stoploss SNV"

msg "Test Completed Successfully"

true

