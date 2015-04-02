#!/bin/sh

#------------------------------------------------------------------------------
# Create a zip file for distribution
# Note: Only binary data is included (no raw gene info / genomes)
#
#                                      Pablo Cingolani 2010
#------------------------------------------------------------------------------

VERSION="3_1"
VERSION_REV=$VERSION""
DIR=$HOME/snpEff_$VERSION_REV
rm -rvf $DIR
mkdir $DIR

# Copy core files
cp snpEff.config snpEff.jar SnpSift.jar $DIR
cp -rvfH galaxy scripts $DIR

cd $DIR
rm -rvf `find . -name "CVS" -type d`
cd -

# Create 'core' zip file
cd $HOME
ZIP="snpEff_v"$VERSION_REV"_core.zip"
rm -f $ZIP 2> /dev/null
zip -r $ZIP snpEff_$VERSION_REV
cd -

# Create ZIP file for each database
for d in `ls data/egran*/snpEffectPredictor.bin`
do
	DIR=`dirname $d`
	GEN=`basename $DIR`
	
	echo $GEN
	ZIP="snpEff_v"$VERSION"_"$GEN".zip"
	zip -r $ZIP data/$GEN/*.bin
done

# Look for missing genomes
echo Missing genomes:
ls -d data/*/snpEffectPredictor.bin | grep -v genomes | cut -f 2 -d / | sort > genomes_bins.txt
ls -d data/* | grep -v genomes | cut -f 2 -d / | sort > genomes_dirs.txt
diff genomes_dirs.txt genomes_bins.txt | grep "^<"

