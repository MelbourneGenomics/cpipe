#!/bin/bash
# vim: ts=4:expandtab:sw=4:cindent
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

function err() {
    echo "======================== Error ========================="
    echo 
    echo  "$1" | fmt -w 60
    echo
    exit 1
}

function usage() {

    echo
    echo "Usage: new_target_region.sh <analysis profile id> <target region bed>"
    echo
    exit 1

}

# Need at least two arguments
[ -z "$2" ] && usage

ANALYSIS_PROFILE_ID="$1"
DIAGNOSTIC_TARGET="$2"

[ -e "$DIAGNOSTIC_TARGET" ] || err "Could not find file $DIAGNOSTIC_TARGET"

[ -e "./bpipe" ] || err "Please run this script from the Cpipe base directory"

[ -e ./designs/flagships/$ANALYSIS_PROFILE_ID.bed ] && \
        err "The analysis profile that you specified already exists in the designs/flagships directory"\
            "To recreate the analysis profile, please first delete it manually from the designs/flagships"\
            "directory"

./bpipe run \
         -p TARGET_REGION_ID="$ANALYSIS_PROFILE_ID" \
         -d designs/$ANALYSIS_PROFILE_ID pipeline/scripts/reformat_bed_files.groovy \
        "$DIAGNOSTIC_TARGET" || err "Creation of target region failed. Please see previous error messages"


[ -e designs/$ANALYSIS_PROFILE_ID/succeeded.txt ] || \
        err "Creation of target region failed. Please see previous error messages"

[ -e designs/$ANALYSIS_PROFILE_ID/$ANALYSIS_PROFILE_ID.bed ] || \
        err "Creation of target region failed. Please see previous error messages"

# Create dummy versions of other required files
cd designs/$ANALYSIS_PROFILE_ID
touch $ANALYSIS_PROFILE_ID.genes.txt
touch $ANALYSIS_PROFILE_ID.transcripts.txt

echo
echo "Succeeded."
echo
