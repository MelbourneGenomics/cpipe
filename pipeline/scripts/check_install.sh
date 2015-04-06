#!/bin/bash
# vim: ts=4:expandtab:sw=4:cindent
############################################################
#
# Installation check script for Melbourne Genomics Pipeline
# dependencies. 
#
# This script does not do the complete job, but it checks
# a number of things and can improve over time.
#
# Author: Simon Sadedin, ssadedin@mcri.edu.au
#
############################################################

# Helper functions
function err() {
    echo
    echo "========================= ERROR =================================="
    echo
    echo "$1" | fmt -w 100
    echo
    echo "=================================================================="
    echo
    exit 1
}

function warn() {
    echo
    echo "================================================================"
    echo "WARNING: $1" | fmt -w 100
    echo "================================================================"
    echo
}


function msg() {
    echo
    echo "================================================================"
    echo "$1"
    echo "================================================================"
    echo
}

#
# Support for "quiet" operation - only prompt if QUIET != false
#
function prompt() {

    PROMPT="$1"
    DEFAULT="$2"
   
    if $QUIET;
    then
        REPLY="$DEFAULT"
    else
        read -p "$1"
    fi
}

function compile() {
    PROGRAM="$1"
    msg "Check $PROGRAM is compiled ..."
    if [ ! -e $PROGRAM ];
    then
            prompt "$PROGRAM does not seem to be compiled: do you want me to compile it? (y/n): " "y"
            if [ "$REPLY" == "y" ];
            then
                pushd `dirname $PROGRAM` || pushd `dirname $(dirname $PROGRAM)`
                [ ! -e Makefile ] && cd ..
                [ ! -e Makefile ] && err "Cannot find Makefile for $PROGRAM"
                make || err "$PROGRAM failed to compile"
                popd
                [ ! -e $PROGRAM ] && err "Could not find $PROGRAM even after compiling it"
            else
                warn "You will need to compile $PROGRAM or edit config.groovy to point to your installation manually"
            fi
    fi
}

QUIET=false
if [ "$1" == "-q" ];
then
  QUIET=true
fi

[ -e pipeline ] || \
        err "I can't see the pipeline directory. Maybe you didn't clone the repository, or you're running this script from the wrong location?"

eval `sed 's/\/\/.*$//' pipeline/config.groovy` 

[ ! -e pipeline/config.groovy ] && \
        err "You haven't created the file pipeline/config.groovy yet: you need to copy the file pipeline/config.groovy.template and edit it."

msg "Check base location is correct ..."
[ ! -e "$BASE/pipeline" ] &&  \
        err "Cannot see $BASE/pipeline - please check BASE is defined correctly in config.groovy. It should probably be "`pwd`

msg "Checking dependencies ..."

msg "Checking version of Java ..."

type $JAVA > /dev/null || \
        err "No Java executable could be located in the current PATH. Please ensure Oracle, Sun or OpenJDK Java is the default Java in your path, and set JAVA_HOME to the corresponding Java location."

$JAVA -version | grep -q gij && \
        err "You are using the GNU Java implementation which is not compatible with the pipeline. Please ensure Oracle, Sun or OpenJDK Java is the default Java in your path, and set JAVA_HOME to the corresponding Java location."

[ -z "$JAVA_HOME" ] && err "The JAVA_HOME environment variable is not defined. Please set it to the location of your Java installation."

compile "$BWA"

compile "$SAMTOOLS/samtools"

compile "$BEDTOOLS/bin/bedtools"

msg "Check GATK is downloaded and available"
[ -e $GATK/GenomeAnalysisTK.jar ] || \
        err "Could not locate $GATK/GenomeAnalysisTK.jar file. Please download and install GATK to $GATK/. See instructions in tools/gatk/README." 

msg "Check Annovar is downloaded and available"
[ -e $ANNOVAR/annotate_variation.pl ] || \
        err "Could not locate Annovar script. Please download and install Annovar to $ANNOVAR/"

msg "Check Annovar database exists"
for i in hg19_snp138.txt hg19_avsift.txt hg19_esp5400_all.txt hg19_refGene.txt hg19_ALL.sites.2010_11.txt hg19_phastConsElements46way.txt; 
do
    [ -e $ANNOVAR/../humandb/$i ] || {
        err "Failed to find all necessary Annovar data files ($i): please use Annovar downdb to download all the data files. See tools/annovar/README for assistance."
    }
done

msg "Check VEP database downloaded for version $VEP_VERSION..."
[ -e $VEP/../vep_cache/homo_sapiens/$VEP_VERSION/1 ] || \
    err "Failed to find downloaded VEP data files. Please install homosapiens cache and FASTA files using: 'cd $VEP; perl INSTALL.pl -c ../vep_cache"

msg "Check that reference FASTA exists"
[ -e "$REF" ] || err "Reference FASTA file $REF could not be found. Please place it there or change config.groovy to point to your reference"

msg "Check reference FASTA is indexed"

[ -e "$REF.fai" ] || err "Reference FASTA file $REF is not indexed. Please run samtools faidx to index it"

[ -e "$REF.bwt" ] || err "Reference FASTA file $REF is not indexed by bwa. Please run 'cd "`dirname $REF`"; ../tools/bwa/0.7.5a/bwa index -a bwtsw ucsc.hg19.fasta' on reference file ($REF) to index it"

[ -e `echo "$REF" | sed 's/\.fa$/.dict/'` ] || err "Reference FASTA file $REF doesn't have a dictionary. Please run Picard CreateSequenceDictionary to make the dictionary (or download the .dict file)."

find `dirname $REF`/ -name '*.bwt' -mtime +180 | grep -q bwt && {
    warn "The BWA index on your reference is more than 180 days old. If you experience errors in the alignment stage, please try re-indexing your data"
    prompt "Press enter to continue" " "
}

[ -e "$DBSNP" ] || err "The DBSNP file $DBSNP does not exist. Please download it."

[ -e "$GOLD_STANDARD_INDELS" ] || {
    if [ -e "$GOLD_STANDARD_INDELS.gz" ];
    then
        err "The indel file $GOLD_STANDARD_INDELS is still gzipped. Please use gunzip to uncompress it."
    else
        err "The indel file $GOLD_STANDARD_INDELS does not exist. Please download it from the GATK resource bundle."
    fi
}

msg "Check ulimit ..."
MAX_OPEN_FILES=`ulimit -n` 
if [ "$MAX_OPEN_FILES" -lt 2048 ]; 
then 
    warn "The limit on open files is set to $MAX_OPEN_FILES. Cpipe may require more open files than this. Consider adding 'ulimit set -S -n 2048' to your .bashrc file."
fi

msg "Success: all the dependencies I know how to check are OK"

