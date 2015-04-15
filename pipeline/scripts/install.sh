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
        read -p "$1 "
    fi
}

#
# Helper function to execute standard compile process for tools that are
# shipped in source form
#
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

function load_config() {
    if [ -z "$BASE" ];
    then
        BASE="."
    fi
    eval `sed 's/\/\/.*$//' $BASE/pipeline/config.groovy` 
}

function set_config_variable() {
    NAME="$1"
    VALUE="$2"
    sed -i 's,'$NAME'=".*$,'$NAME'="'$VALUE'",g' $BASE/pipeline/config.groovy
    load_config
}

echo '
========================================

 #####   ######   ###  ######   #######  
#     #  #     #   #   #     #  #        
#        #     #   #   #     #  #        
#        ######    #   ######   #####    
#        #         #   #        #        
#     #  #         #   #        #        
 #####   #        ###  #        #######  

========================================
'

QUIET=false
if [ "$1" == "-q" ];
then
  QUIET=true
fi

[ -e pipeline ] || \
        err "I can't see the pipeline directory. Maybe you didn't clone the repository, or you're running this script from the wrong location?"

[ ! -e pipeline/config.groovy ] && {
        prompt "You haven't created the file pipeline/config.groovy yet. Do you want me to copy this file from the default template for you? (y/n)" "y"
        if [ "$REPLY" == "y" ]
        then
            cp -uv pipeline/config.groovy.template pipeline/config.groovy
        fi

        BASE=`pwd`
        set_config_variable BASE "$BASE"
        load_config
}


load_config

msg "Check base location is correct ..."
[ ! -e "$BASE/pipeline" ] &&  \
        err "Cannot see $BASE/pipeline - please check BASE is defined correctly in config.groovy. It should probably be: "`pwd`

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
[ -e $GATK/GenomeAnalysisTK.jar ] || {
    echo "
 The recommended version of GATK for Cpipe is 3.3. However license 
 terms prevent Cpipe from including this version with Cpipe. Cpipe
 includes GATK 2.3.9 which can be used instead. If you wish to use
 a later version of GATK, please abort this script (Ctrl-c), download
 that version after separately agreeing to the license terms, and 
 place the jar file in tools/gatk/<version>/GenomeAnalysisTK.jar. Once
 you have done that, set the GATK variable in pipeline/config.groovy 
 appropriately and re-run this script.
    "
    prompt "Continue with GATK 2.3.9? (y/n)" "y"
    if [ "$REPLY" == "y" ];
    then
        set_config_variable GATK "$BASE/tools/gatk/2.3.9"
        set_config_variable GATK_LEGACY "true"
    else
        msg "WARNING: your installation will not work unless you set GATK manually youself in pipeline/config.groovy"
    fi
}

msg "Check Annovar is downloaded and available"
[ -e $ANNOVAR/annotate_variation.pl ] || {

  pushd $BASE/tools/annovar > /dev/null

  ANNOVAR_DOWNLOAD=`ls -t *.tar.gz | head -1`

  if [ ! -e "$ANNOVAR_DOWNLOAD" ];
  then
        echo "
  Due to license restrictions, Annovar cannot be included
  with Cpipe. If you wish to use Annovar, please download it
  separately and place the downloaded file (tar.gz) in the 
  following location: $BASE/tools/annovar

  Once you have done this, press enter to continue.
      "
      read  -p "Press enter when you have placed the Annovar tar.gz file in $BASE/tools/annovar ..."
  fi

  ANNOVAR_DOWNLOAD=`ls -t *.tar.gz | tail -1`
  if [ -e "$ANNOVAR_DOWNLOAD" ];
  then
      tar -xzf $ANNOVAR_DOWNLOAD
      ANNOVAR_PL=`find . -iname annotate_variation.pl | xargs ls -t | head -1` 
      ANNOVAR_VERSION=`perl $ANNOVAR_PL | grep Version | grep -o '[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}'`
      ANNOVAR_DIR=`dirname $ANNOVAR_PL`
      mv $ANNOVAR_DIR $ANNOVAR_VERSION

      set_config_variable ANNOVAR "$BASE/tools/annovar/$ANNOVAR_VERSION"
  else
      err "Could not find tar.gz file in $BASE/tools/annovar"
  fi
  popd > /dev/null
}

msg "Check Annovar database exists"
MISSING_ANNOVAR=""
ANNOVAR_DB_FILES="hg19_snp138.txt hg19_avsift.txt hg19_esp5400_all.txt hg19_refGene.txt hg19_ALL.sites.2010_11.txt hg19_phastConsElements46way.txt"
for i in $ANNOVAR_DB_FILES ;  
do
    [ -e $ANNOVAR/../humandb/$i ] || {
        MISSING_ANNOVAR="$i $MISSING_ANNOVAR"
    }
done

if [ -z "$MISSING_ANNOVAR" ];
then
    echo "
    One or more Annovar database files is not present. Do you want to 
    download Annovar files now using the built in download script?
    NOTE: downloading may take some time and cause large files to be
    downloaded. Please try to avoid having this process be interrupted.

    Missing files:

    $MISSING_ANNOVAR
    "

    prompt "Download Annovar files? (y/n)" "y"

    if [ "$REPLY" == "y" ];
    then
        $BASE/pipeline/scripts/download_annovar_db.sh $ANNOVAR $BASE/tools/annovar/humandb \
            || err "Failed to download Annovar databases"
    else
        msg "WARNING: Cpipe will not operate correctly if Annovar database file are not present."
    fi
fi

msg "Check VEP database downloaded for version $VEP_VERSION..."
[ -e $VEP/../vep_cache/homo_sapiens/$VEP_VERSION/1 ] || {
    echo "
    Cpipe uses the Variant Effect Predictor from Ensembl
    to perform annotation of variants.
    
    To work, VEP must first be installed, and then reference
    files must be downloaded. Would you like to launch the
    VEP installer script now?

    Please answer "y" when asked to install cache files, and choose
    either homo_sapiens_vep_74.tar.gz or homo_sapiens_refseq_vep_74.tar.gz
    to download. Note: Cpipe is tested and developed using homo_sapiens_vep_74.tar.gz.
    "

    prompt "Do you want to run the VEP installer now? (y/n)" "y"
    if [ "$REPLY" == "y" ];
    then
        cd $VEP; 
        perl INSTALL.pl -c ../vep_cache || err "Failed to run VEP installer"
    else
        msg "WARNING: Cpipe will not operate correctly if VEP is not installed"
    fi
}

msg "Configuring Condel Plugin ..."

sed -i 's,do not use,'$CONDEL/config',' $CONDEL/config/condel_SP.conf || err "Unable to configure Condel plugin"

msg "Check that reference FASTA exists"
[ -e "$REF" ] ||  {
  echo "
  It appears that you have not set or downloaded a human
  genome reference yet ($REF is not found). An appropriate reference
  and associated files can be downloaded from the GATK reference bundle.
  
  "
  prompt "Would you like to download and install files from the GATK bundle now? (y/n)" "y"

  GATK_BUNDLE_FILES=`basename $DBSNP`.gz" "`basename $GOLD_STANDARD_INDELS`.gz" ucsc.hg19.dict.gz ucsc.hg19.fasta.gz ucsc.hg19.fasta.fai.gz"
  if [ "$REPLY" == "y" ];
  then
      pushd $REFBASE /dev/null
      for f in $GATK_BUNDLE_FILES ;
      do
          if [ ! -e `echo $f | sed 's/.gz$//'` ];
          then
              BUNDLE_URL="ftp://ftp.broadinstitute.org/bundle/2.8/hg19/$f" 
              wget --user=gsapubftp-anonymous \
                   --password=cpipe.user@cpipeline.org \
                   $BUNDLE_URL || err "Unable to download file $BUNDLE_URL"

              gunzip $f || err "Failed to unzip downloaded file $f"
          else
             echo "File $f already exists ..."
          fi
      done 
      popd > /dev/null
  else
      msg "WARNING: Cpipe will not operate correctly if reference files are not available"
  fi
}

msg "Check reference FASTA is indexed"

[ -e "$REF.fai" ] || err "Reference FASTA file $REF is not indexed. Please run samtools faidx to index it"

[ -e "$REF.bwt" ] || {

    prompt "Reference FASTA file $REF is not indexed by bwa. Do you want to index it now? (y/n)?" "y"

    if [ "$REPLY" == "y" ];
    then
        cd "`dirname $REF`"; 
        $BWA index -a bwtsw `basename $REF` || err "Indexing reference $REF using bwa failed"
    else
        msg "WARNING: Cpipe will not operate correctly unless the reference genome is indexed"
    fi
}

[ -e `echo "$REF" | sed 's/\.fa$/.dict/'` ] \
        || err "Reference FASTA file $REF doesn't have a dictionary. Please run Picard CreateSequenceDictionary to make the dictionary (or download the .dict file)."


find `dirname $REF`/ -name '*.bwt' -mtime +180 | grep -q bwt && {
    warn "The BWA index on your reference is more than 180 days old. If you experience errors in the alignment stage, please try re-indexing your data"
    prompt "Press enter to continue" " "
}

msg "Check if genome file $HG19_CHROM_INFO exists ..."
[ -e "$HG19_CHROM_INFO" ] || {
    prompt "The chromosome size file (.genome) does not exist. Create it by downloading from UCSC (requires MySQL client) (y/n)" "y"
    if [ "$REPLY" == "y" ];
    then
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
            "select chrom, size from hg19.chromInfo"  > $HG19_CHROM_INFO || err "Failed to download genome file"
    else
        msg "WARNING: Cpipe will not operate correctly unless the $HG19_CHROM_INFO file is created"
    fi
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

