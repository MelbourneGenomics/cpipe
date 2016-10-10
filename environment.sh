#!/usr/bin/env bash

# Util functions
function join() {
    local IFS=$1
    shift
    echo "$*"
}


# Work out the directory name
ROOT=$(readlink -f $(dirname $BASH_SOURCE))

# Load config groovy
source $ROOT/pipeline/scripts/config_groovy_util.sh
load_config

# Load virtualenv
source ${TOOLS}/python/bin/activate

# Add all tool directories and bin folders to PATH
export PATH=`join ':' $TOOLS/*/`:`join ':' $TOOLS/*/bin/`:$PATH
export HTSLIB_DIR=$TOOLS/htslib
export PERL5LIB=$TOOLS/perl_lib/lib/perl5:$TOOLS/perl/lib:$TOOLS/vep_plugins:${TOOLS}/vep_libs
export C_INCLUDE_PATH=${C_INCLUDE_PATH}:`echo ${TOOLS}/c_libs/*`
export PERL_SRC=${TOOLS}/perl
