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
#     Runs all unit tests 
#
########################################################

# Load paths from config.groovy
source $PWD/pipeline/scripts/config_groovy_util.sh
load_config

# Variables
GROOVY_TEST_LIBRARIES="$PWD/tools/groovy/lib/*:$PWD/tools/java_libs/*"
GROOVYC="${GROOVY}c"

# Python tests
pushd pipeline/scripts
    python -m unittest discover -s ../tests -p '*_test.py' -v
popd

# Groovy tests
pushd pipeline/tests
    # compile scripts and tests
    mkdir -p tmp
    sh $GROOVYC -cp $GROOVY_TEST_LIBRARIES -d tmp ../scripts/*.groovy ./*.groovy
    $JAVA -cp "$GROOVY_TEST_LIBRARIES:$PWD/tmp" RunAll
popd
