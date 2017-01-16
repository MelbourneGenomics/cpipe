#!/usr/bin/env bash

# Util functions
function join() {
    local IFS=$1
    shift
    echo "$*"
}

# Work out the directory name
if [ -n "$ZSH_VERSION" ]; then
    ROOT=$(readlink -f $(dirname $0))
else
    ROOT=$(readlink -f $(dirname $BASH_SOURCE))
fi

# Add all tool directories and bin folders to PATH
export SYS_JAVA=`which java` # Export the old system java before we override it
export PATH=${TOOLS}/bin:${TOOLS}/maven/bin:${TOOLS}/bpipe/bin:${PATH}
export HTSLIB_DIR=$TOOLS/lib
export PERL5LIB=$TOOLS/perl_lib/lib/perl5:$TOOLS/perl/lib:$TOOLS/vep_plugins:${TOOLS}/vep_libs
export CPATH="${TOOLS}/include"
export CFLAGS="$CFLAGS -I${TOOLS}/include"
export CPPFLAGS="$CPPFLAGS -I${TOOLS}/include"
export LDFLAGS="$LDFLAGS -L${TOOLS}/lib"
export LD_LIBRARY_PATH=${TOOLS}/lib:${LD_LIBRARY_PATH}
export JAVA_OPTS #Pass JAVA_OPTS to the script in c_libs/bin/java
export TMPDIR #TMPDIR is set in config.groovy.
export CPIPE_ROOT=${ROOT}

# Load config groovy
source ${ROOT}/lib/config_groovy_util.sh
load_config

# Load virtualenv
source ${TOOLS}/python/bin/activate

