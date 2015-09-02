#!/bin/bash
# vim: expandtab:shiftwidth=4:ts=4:cindent
###########################################################################
#
# This file is part of Cpipe.
# 
# Cpipe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, under version 3 of the License, subject
# to additional terms compatible with the GNU General Public License version 3,
# specified in the LICENSE file that is part of the Cpipe distribution.
#
# Cpipe is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.  If not, see <http:#www.gnu.org/licenses/>.
#
###########################################################################

function err() {
    echo
    echo "========================== ERROR =============================="
    echo "$1" | fmt -w 80
    echo "==============================================================="
    usage
    exit 1
}

function abs_path() {
    echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
}

function usage() {
    echo '
Usage: create_batch.sh <batch identifier> <target region identifier> <exome bed file>

1. The batch identifier should be an existing directory that you have 
   created under 'batches', already containing a data subdirectory containign 
   FASTQ files.

2. The target region identifier should be one of the predefined sub target regions 
   that is set up in the designs directory. These are the region for
   analysis within the larger exome.

3. The exome BED file defines the regions covered by the whole exome.

Example:

    ./pipeline/scripts/create_batch.sh cardiac_001 CARDIAC101 designs/CARDIAC101/CARDIAC101.bed

Note in this example, the region for analysis is the same as the whole design region
captured by the kit.

'
}

echo '
=================== Cpipe Batch Creation Helper Script ========================
'


[ ! -e batches ] && 
    err "Please run this script from the root directory of the pipeline distribution (the parent directory of pipeline, batches, etc)"

[ -z "$1" ] && err "Please specify the batch identifier"

BATCH_ID="$1"

# Allow specification of batch by its directory
[ ! -e "./batches/$BATCH_ID" ] &&  {
    if [ -e $BATCH_ID ] && [ -e batches/`basename "$BATCH_ID"` ];
    then
        BATCH_ID=`basename $BATCH_ID`
    fi
}


[ ! -e "./batches/$BATCH_ID" ] && 
    err "Could not see a directory called ./batches/$BATCH_ID - have you created the batch directory yet?"

[ -z "$2" ] && 
    err "Please specify the target region identifier, eg: CARDIAC101"

TARGET="$2"
[ -e "designs/${TARGET}/${TARGET}.bed" ] || {
    if [ -e $TARGET ] && [ -e designs/`basename $TARGET` ];
    then
        TARGET=`basename $TARGET`
    fi
}

[ -e "designs/${TARGET}/${TARGET}.bed" ] ||
    err "No file called ${TARGET}.bed could be found in designs/${TARGET}. Is this target region configured yet?"

eval `sed 's/\/\/.*$//' pipeline/config.groovy`


# Default: exome target is the same as diagnostic target / analysis profile
EXOME_TARGET="$3"
if [ -z "$EXOME_TARGET" ];
then
    EXOME_TARGET=$TARGET.bed
fi

# The exome target can be a BED file inside the analysis profile directory
if [ -e "$BASE/designs/${TARGET}/$EXOME_TARGET" ];
then
    EXOME_TARGET="$BASE/designs/${TARGET}/$EXOME_TARGET"
elif [ "./designs/${TARGET}/"`basename $EXOME_TARGET`=="./designs/${TARGET}/$EXOME_TARGET" ] && [ -e "./designs/${TARGET}/$EXOME_TARGET" ];
then
    EXOME_TARGET="$BASE/designs/${TARGET}/"`basename $EXOME_TARGET`
fi

[ -f "$EXOME_TARGET" ] ||
    err "No file called ${EXOME_TARGET} could be found. You may need to specify this file with an absolute path"

# Ensure EXOME_TARGET is expressed as an absolute path
EXOME_TARGET=`abs_path "$EXOME_TARGET"`

cd "batches/$BATCH_ID"

echo '
EXOME_TARGET="'$EXOME_TARGET'"
' > target_regions.txt

[ -e samples.txt ] && err \
        "A samples.txt file already exists for this batch. Please remove it and any previous analysis manually before rerunning this script."

# echo "$GROOVY -cp $BASE/tools/groovy-ngs-utils/1.0/groovy-ngs-utils.jar $BASE/pipeline/scripts/files_to_sample_info.groovy -batch $BATCH_ID -disease $TARGET data/*.fastq.gz > samples.txt"

echo "Creating sample meta data file batches/$BATCH_ID/samples.txt with disease $TARGET..."

$GROOVY -cp $BASE/tools/groovy-ngs-utils/1.0/groovy-ngs-utils.jar $BASE/pipeline/scripts/files_to_sample_info.groovy -batch $BATCH_ID -disease $TARGET data/*.fastq.gz > samples.txt || \
        err "Configuring the samples.txt file failed. Please check error messages for the reason and ensure your FASTQ files are named correctly"


mkdir analysis || err "Unable to create analysis directory"

echo "

Done.

To run this analysis, change directory to ./batches/$BATCH_ID/analysis and execute:

../../../bpipe run ../../../pipeline/pipeline.groovy ../samples.txt

"
