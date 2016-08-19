#!/usr/bin/env bash

ROOT=`readlink -f $(dirname $0)`

# Load config variables
source $ROOT/pipeline/scripts/load_config_groovy.sh

function join() {
    local IFS=$1
    shift
    echo "$*"
}

function compile {
# $1 is the directory to compile in
    pushd $1

    if [[ -f  configure.ac ]]; then
        autoconf
    fi
    if [[ -f configure ]]; then
        yes | ./configure
    fi
    if [[ -f Configure ]]; then
        yes | ./Configure -d
    fi
    if [[ -f Makefile ]]; then
        make
    fi

    popd
}

# Compile everything
for DIRECTORY in $TOOLS_ROOT/*/; do
    compile $DIRECTORY
done

# Add all tool directories and bin folders to PATH
export PATH=`join ':' $TOOLS_ROOT/*/`:`join ':' $TOOLS_ROOT/*/bin/`:$PATH
export HTSLIB_DIR=$TOOLS_ROOT/htslib

# Get HTSlib module
cpanm Bio::DB::HTS