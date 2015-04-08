#!/bin/bash
# vim: ts=4:expandtab:sw=4
########################################################
#
#      Melbourne Genomics Health Alliance
#
########################################################

EMAILS="ssadedin@gmail.com"

TESTS='*_test'

if [ ! -z "$1" ];
then
    TESTS="$1"
fi

#
# Helper functions
function err() {
        (
        echo
        echo "========================= ERROR =================================="
        echo
        echo "$1" | fmt -w 80
        echo
        echo "=================================================================="
        echo
        ) | tee error.log

        mail -s "WARNING: Melbourne Genomics SelfTest Failure" $EMAILS  < error.log

        exit 1
}

function msg() {
        echo
        echo "================================================================"
        echo "$1" | fmt -w 80
        echo "================================================================"
        echo
}

function annovar_has() {
  # TODO: make this more robust with a real CSV parser, etc
  [ -z "$ANNOVAR_CSV" ] && err "ANNOVAR_CSV variable not defined"
  CHR=$1
  POS=$2
  TYPE="$3"
  printf "Check $TYPE at $CHR:$POS "

  cat $ANNOVAR_CSV | grep '"'$CHR'"' | grep '"'$POS'"' | grep -q '"'"$TYPE"'"' \
    || err "Annovar file $ANNOVAR_CSV did not have expected variant of type '$TYPE' at  $CHR:$POS"

  echo "PASS" 
}

# This script should ALWAYS be launched from the main pipeline root directory
if [ ! -e pipeline ] && [ ! -e designs ];
then
    err "Could not find pipeline directory. Please run this script from the root of the pipeline distribution."
fi

# Load all the pipeline configuration settings
eval `sed 's/\/\/.*$//' pipeline/config.groovy` 

# Now run each test
for i in pipeline/tests/$TESTS ;
do
    pushd $i > /dev/null
    
    #echo "Would run $i"
    source run.sh

    popd > /dev/null
done

