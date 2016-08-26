#!/usr/bin/env bash

### Logging ####
CURRENT_DIR=$(readlink -f $(dirname $BASH_SOURCE))
LOG_FILE=$CURRENT_DIR/install.log
echo "## Running installer on `date` ##" > $LOG_FILE
echo "## Logging to $LOG_FILE ##"

source $CURRENT_DIR/common/compile.sh #Load compile function
ROOT=$CURRENT_DIR/..

# Load config variables
source $ROOT/pipeline/scripts/config_groovy_util.sh
set_config_variable BASE $(readlink -f $ROOT)

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

#Fastqc is a perl script so needs to be chmod +x'd
chmod +x $TOOLS/fastqc/fastqc

# Get HTSlib module
cpanm Bio::DB::HTS
