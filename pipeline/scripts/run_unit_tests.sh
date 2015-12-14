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

export GROOVY="$PWD/tools/groovy/2.3.4/bin/groovy"
export GROOVY_NGS="$PWD/tools/groovy-ngs-utils/1.0.2"
export EXCEL="$PWD/tools/excel/1.0"
GROOVY_TEST_LIBRARIES="$PWD/tools/groovy/2.3.4/lib/groovy-2.3.4.jar:$PWD/tools/groovy/2.3.4/lib/hamcrest-core-1.3.jar:$PWD/tools/groovy/2.3.4/lib/junit-4.11.jar:$PWD/pipeline/tests/lib/cpsuite-1.2.6.jar:$EXCEL/excel.jar:$GROOVY_NGS/groovy-ngs-utils.jar"
GROOVYC="$PWD/tools/groovy/2.3.4/bin/groovyc"

# python tests
pushd pipeline/tests
python -m unittest discover -s . -p '*_test.py' -v
popd

# groovy tests
pushd pipeline/tests
# compile scripts and tests
mkdir -p tmp
sh $GROOVYC -cp $GROOVY_TEST_LIBRARIES -d tmp ../scripts/*.groovy ./*.groovy
java -cp $GROOVY_TEST_LIBRARIES:tmp RunAll
popd
