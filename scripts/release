#!/bin/bash
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
# Release Helper Script - This script is designed to assist in versioning
# Cpipe. 
#
###########################################################################

function err() {
    echo 
    echo "ERROR: $1"
    echo
    exit 1
}

if [ ! -e pipeline/scripts/release.sh ];
then
    err "Please run this script from the root of the pipeline distribution."
fi

echo
echo "================================================="
echo
echo " Cpipe Release Script"
echo
echo "================================================="
echo "
What this script does:

  1  allows you to enter a version number
  2  allows you to enter release notes for the release, based
     on changelog from git
  3  creates a tag for the version in the git repository
  4  pushes tag to main repository

How it works:

  The version number is stored in a file in the top level directory:

    version.txt

  This file contains the current version of the pipeline, and this is
  what is imprinted into the provenance files that are output
  when the pipeline is run. This script both updates the version number
  and creates a tag in git reflecting the version, and pushes that
  to the main github repository.

"

echo "The previous released version was: "
echo
cat version.txt
echo

read -p "Please enter the version number for the release: "

VERSION=$REPLY
echo
echo "Version is $VERSION"
echo
echo "
    Press enter now to open release notes editor. The editor will open with the 
    git changelog and past release notes in separate panes of the editor.
"

read

PREV_VERSION=`cat version.txt`

echo "==================================================

Melbourne Genomics Demonstration Project Pipeline

Version $VERSION Release Notes

==================================================

# Please enter your release notes here. The pane
# opposite has a list of changes committed since the 
# last release was made. Please delete this comment.

" > notes.tmp.txt

git log ${PREV_VERSION}..HEAD | grep -v '^Author' | grep -v '^commit' | grep -v '^Date' >> git.tmp.txt

touch .notes.tmp.txt.timestamp

vim -O notes.tmp.txt git.tmp.txt 

if [ .notes.tmp.txt.timestamp -nt notes.tmp.txt ];
then
        echo "Release notes not saved - aborting ...."
        echo
        exit 1
fi

cp notes.tmp.txt docs/ReleaseNotes-$VERSION.txt || err "Unable to copy release notes to docs directory"

git add docs/ReleaseNotes-$VERSION.txt || err "Unable to add new release notes to git repo"
echo $VERSION > version.txt || err "Unable to update version file"
git commit -m "Release notes for version $VERSION" version.txt docs/ReleaseNotes-$VERSION.txt || err "Unable to commit to git repository"
git tag $VERSION || err "Unable to create tag for version $VERSION"
git push || err "Unable to push code to origin repo"
git push origin $VERSION || err "Unable to push tag $VERSION to origin repo"

rm .notes.tmp.txt.timestamp notes.tmp.txt 

echo
echo "Done."
echo

