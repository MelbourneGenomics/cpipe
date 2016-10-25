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

CROOT=${TOOLS}/c_libs

# Load virtualenv
source ${TOOLS}/python/bin/activate

# Add all tool directories and bin folders to PATH
export SYS_JAVA=`which $JAVA` # Export the old system java before we override it 
export PATH=${TOOLS}/bin:${TOOLS}/maven/bin:${TOOLS}/bpipe/bin:${PATH}
export HTSLIB_DIR=$TOOLS/htslib
export PERL5LIB=$TOOLS/perl_lib/lib/perl5:$TOOLS/perl/lib:$TOOLS/vep_plugins:${TOOLS}/vep_libs
export CPATH="${TOOLS}/c_libs/include"
export CFLAGS="$CFLAGS -I${TOOLS}/c_libs/include"
export CPPFLAGS="$CPPFLAGS -I${TOOLS}/c_libs/include"
export LDFLAGS="$LDFLAGS -L${TOOLS}/c_libs/lib"
export LD_LIBRARY_PATH=${TOOLS}/c_libs/lib:${LD_LIBRARY_PATH}
export JAVA_OPTS #Pass JAVA_OPTS to the script in c_libs/bin/java
