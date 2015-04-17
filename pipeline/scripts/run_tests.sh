#!/bin/bash
# vim: ts=4:expandtab:sw=4
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
#
#      Selftest Script
#
########################################################

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

        if [ ! -z "$EMAILS" ];
        then
            mail -s "WARNING: Melbourne Genomics SelfTest Failure" $EMAILS  < error.log
        fi

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
    source ./run.sh

    popd > /dev/null
done

