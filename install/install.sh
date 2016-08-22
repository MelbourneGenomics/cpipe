#!/usr/bin/env bash

### Logging ####
LOG_FILE=$(readlink -f $(dirname $0)/install.log)
echo "## Running installer on `date` ##" > $LOG_FILE
echo "## Logging to $LOG_FILE ##"

source `dirname $0`/common/compile.sh
ROOT=`readlink -f $(dirname $0)/..`
TOOLS_ROOT=`readlink -f $(dirname $0)/..`

# Load config variables
source $ROOT/pipeline/scripts/load_config_groovy.sh

function join() {
    local IFS=$1
    shift
    echo "$*"
}

# Compile everything
for DIRECTORY in $TOOLS/*/; do
    echo -n "Compiling `basename $DIRECTORY`..."
    compile $DIRECTORY &>> $LOG_FILE
    if [[ $? -ne 0 ]] ; then
        echo "Not required."
    else
        echo "Done."
    fi
done

# Add all tool directories and bin folders to PATH
export PATH=`join ':' $TOOLS_ROOT/*/`:`join ':' $TOOLS_ROOT/*/bin/`:$PATH
export HTSLIB_DIR=$TOOLS_ROOT/htslib

# Get HTSlib module
cpanm Bio::DB::HTS