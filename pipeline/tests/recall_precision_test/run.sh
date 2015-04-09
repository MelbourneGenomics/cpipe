#!/bin/bash
# vim: ts=4:expandtab:sw=4:cindent
###########################################################
#
# Build a batch directory containing our test data 
#

RUN=true
if [ "$1" == "norun" ];
then
    RUN=false
    source go.sh
fi

HERE=`pwd`

# The variant database we will use - we do not want to interfere with any 
# previous variant database
VARIANT_DB=$HERE/variants.db

BATCH_DIR=$BASE/batches/recall_precision_test

cd $BASE

if $RUN ;
then
    if [ -e $BATCH_DIR ];
    then
        rm -rf $BATCH_DIR
    fi


    cd $HERE
    [ -z "$EXOME_TARGET" ] && err "This test requires an EXOME_TARGET to be defined in the main config.groovy file"

    if [ -z "data/NA12878CHR22_L001_R1.fastq.gz" ] || [ -z "data/NA12878CHR22_L001_R2.fastq.gz" ] ;
    then 
        err "The test data required for the precision / recall test is not available. Please download it."
    fi

    if [ ! -e NA12878.target.vcf ];
    then 
        msg "Fetching gold standard calls ..."

        if [ ! -e data/NA12878.vcf.gz ];
        then
            err "Please download a set of gold standard hg19 NA12878 calls and place them in $HERE/data/NA12878.vcf.gz"
        fi
        
        gunzip -c data/NA12878.vcf.gz > data/NA12878.vcf

        msg "Filtering gold standard regions ..."
        $JAVA -Xmx2g -jar $GATK/GenomeAnalysisTK.jar \
             -R $REF \
             -T SelectVariants \
             --variant data/NA12878.vcf \
                 -sn NA12878 \
                 -ef \
             -L $EXOME_TARGET \
             -L chr22 \
             -isr INTERSECTION \
             -o NA12878.target.vcf || err "Unable to filter gold standard variant calls"
    fi
    cd $BASE

    msg "Creating directories ..."
    mkdir -p $BATCH_DIR/data || err "Unable to create batch directory"

    msg "Copying test data ..."
    cp -v $HERE/data/*.fastq.gz $BATCH_DIR/data || err "Unable to copy data to batch directory"

    # Create an appropriate design for our cut down test from the exome target regions
    EXOME_DESIGN=`basename $EXOME_TARGET | sed 's/.bed//g'` 
    [ -e $BASE/designs/${EXOME_DESIGN}_CHR22 ] || {
        grep '^chr22' $EXOME_TARGET > $HERE/data/${EXOME_DESIGN}_CHR22.bed
        ./pipeline/scripts/new_target_region.sh ${EXOME_DESIGN}_CHR22 $HERE/data/${EXOME_DESIGN}_CHR22.bed \
            || err "Unable to create target region for chromsome 22 test data set"
    }

    #[ ! -e designs/$EXOME_DESIGN/$EXOME_DESIGN.bed ] && \
    #    err "This test requires a design be created for the whole exome, but we could not see designs/$EXOME_DESIGN/$EXOME_DESIGN.bed please use pipeline/scripts/new_target_region.sh to create it"

    msg "Creating batch using $EXOME_TARGET analysis profile ..."
    ./pipeline/scripts/create_batch.sh `basename $BATCH_DIR` ${EXOME_DESIGN}_CHR22 \
        || err "Failed to create batch $BATCH_DIR with design ${EXOME_DESIGN}_CHR22"
fi

pushd $BATCH_DIR/analysis

if $RUN ;
then
    msg "Running pipeline ... "
    ../../../bpipe run -p CHECK_FASTQC_FAILURES=false -p VARIANT_DB=$VARIANT_DB ../../../pipeline/pipeline.groovy ../samples.txt | tee output.log
fi

popd

cd "$HERE"

msg "Comparing VCF files  NA12878.target.vcf  to $BATCH_DIR/analysis/variants/NA*.recal.vcf"

# Run VCFSimilarity - this only compares SNVs - we'll do better one day
$JAVA -Xmx4g -cp $GROOVY_HOME/embeddable/groovy-all-2.3.4.jar:$GROOVY_NGS/groovy-ngs-utils.jar VCFSimilarity \
    -s shared.tsv NA12878.target.vcf $BATCH_DIR/analysis/variants/NA*.recal.vcf > similarity.txt \
    || err "Failed to run VCF similarity comparison"

# Check the output - we expect the difference between the samples to be 
$GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar -e \
    'x = Matrix.load("shared.tsv"); println(Math.round(100*(1 - (Math.abs(x[0][0] - x[0][1]) / x[0][1]))))' > sensitivity.txt

SENSITIVITY=`cat sensitivity.txt` 
MIN_SENSITIVITY=90

[  "$SENSITIVITY" -gt $MIN_SENSITIVITY ] || err "Sensitivity = $SENSITIVITY < $MIN_SENSITIVITY"

msg "Test succeeded: sensitivity observed = $SENSITIVITY"





