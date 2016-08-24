#!/usr/bin/env bash

### Logging ####
CURRENT_DIR=$(readlink -f $(dirname $BASH_SOURCE))
LOG_FILE=$CURRENT_DIR/install.log
echo "## Running installer on `date` ##" > $LOG_FILE
echo "## Logging to $LOG_FILE ##"

source $CURRENT_DIR/common/compile.sh #Load compile function
ROOT=$CURRENT_DIR/..

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

sudo chmod +x $TOOLS/fastqc/fastqc

# Add all tool directories and bin folders to PATH
export PATH=`join ':' $TOOLS/*/`:`join ':' $TOOLS/*/bin/`:$PATH
export HTSLIB_DIR=$TOOLS/htslib
export PERL5LIB=$TOOLS/perl_lib/lib/perl5

# Get HTSlib module
cpanm Bio::DB::HTS
