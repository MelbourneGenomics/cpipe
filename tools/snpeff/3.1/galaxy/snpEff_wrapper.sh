#!/bin/sh
#
# Wrapper script for snpEff galaxy tool
#
# Usage: snpEff_wrapper OPTIONS -stats <stats> --output <summary> --genesFile <genes>
#
# OPTIONS (other than those below) are passed directly to snpEFF
#
# -stats NAME:      specifies name of output stats file
# --output NAME:    specifies name of output HTML summary file
# --genesFile NAME: specifies name of output genes file
#
# Edit SNPEFF_JAR and SNPEFF_CONFIG for your local setup
SNPEFF_JAR=/scratch/galaxy/galaxy_dist/tools/snpEff/snpEff.jar
SNPEFF_CONFIG=/scratch/galaxy/galaxy_dist/tools/snpEff/snpEff.config
#
# Process command line
input=
output=
snpEff_args=
while [ ! -z "$1" ] ; do
    case $1 in
	--input)
	    # Name of input file
	    shift; input=$1
	    ;;
	--output)
	    # Name of HTML summary output file (i.e. stdout)
	    shift; output=$1
	    ;;
	-stats)
	    # Name of output stats summary file
	    # This over-rides the -stats option for snpEff
	    shift; statsFile=$1
	    ;;
	--genesFile)
	    # Name of output genes.txt file
	    shift; genesFile=$1
	    ;;
	*)
	    # Collect any other arguments to pass
	    # directly to snpEff
	    snpEff_args="$snpEff_args $1"
	    ;;
    esac
    # Move to next argument
    shift
done
#
# Local versions for stats and genes file
# At the end these will be moved to the names supplied by Galaxy
localStatsFile=snpEff_stats_file
localGenesFile=$localStatsFile.genes.txt
#
# Run snpEff
#
# You can change the amount of memory used by snpEff, just change the -Xmx parameter
# (e.g. use -Xmx2G for 2Gb of memory)
snpEff_cmd="java -Xmx6G -jar $SNPEFF_JAR -c $SNPEFF_CONFIG $snpEff_args -stats $localStatsFile $input"
echo $snpEff_cmd
$snpEff_cmd > $output
#
# Capture the stats and genes.txt files
if [ -f $localStatsFile ] ; then
    echo Moving $localStatsFile to $statsFile
    /bin/mv $localStatsFile $statsFile 2>&1
else
    echo ERROR no stats file found
fi
if [ -f $localGenesFile ] ; then
    echo Moving $localGenesFile to $genesFile
    /bin/mv $localGenesFile $genesFile 2>&1
else
    echo ERROR no genes file found
fi
##
#
