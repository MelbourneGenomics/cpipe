###########################################################
#
# Build a batch directory containing our test data
#
HERE=`pwd`

# The variant database we will use - we do not want to interfere with any
# previous variant database
VARIANT_DB=$HERE/variants.db
BATCH_DIR=$BASE/batches/trio_test

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

# add trio details
python ./pipeline/scripts/update_metadata.py --sample_id NA12877 --name Sex --value "Female" --target $BATCH_DIR/samples.txt
python ./pipeline/scripts/update_metadata.py --sample_id NA12878 --name Sex --value "Male" --target $BATCH_DIR/samples.txt
python ./pipeline/scripts/update_metadata.py --sample_id NA12879 --name Pedigree_File --value "fid000=NA12877,NA12878" --target $BATCH_DIR/samples.txt

pushd $BATCH_DIR/analysis

msg "Running trio pipeline..."

../../../bpipe run -p CHECK_FASTQC_FAILURES=false -p VARIANT_DB=$VARIANT_DB ../../../pipeline/pipeline.groovy ../samples.txt | tee output.log

msg "Checking results ..."

msg "Test Completed Successfully"

popd

true
