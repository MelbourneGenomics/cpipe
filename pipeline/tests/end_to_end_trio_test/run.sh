###########################################################
#
# Build a batch directory containing our test data
#
HERE=$(dirname ${BASH_SOURCE})

BATCH_DIR=$BASE/batches/end_to_end

cd $BASE
if [ -e $BATCH_DIR ];
then
    rm -rf $BATCH_DIR
fi

msg "Creating directories ..."
mkdir -p $BATCH_DIR/data || err "Unable to create batch directory"

msg "Copying test data ..."
cp -v $HERE/data/*.fastq.gz $BATCH_DIR/data || err "Unable to copy data to batch directory"

msg "Creating batch using ALL analysis profile ..."
./pipeline/scripts/create_batch.sh `basename $BATCH_DIR` "ALL" || err "Failed to create batch"

pushd $BATCH_DIR/analysis

msg "Running pipeline ... "

../../../bpipe run -p CHECK_FASTQC_FAILURES=false -p VARIANT_DB=$VARIANT_DB ../../../pipeline/pipeline.groovy ../samples.txt | tee output.log

msg "Checking results ..."


