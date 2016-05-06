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

CUSTOM_ASSETS_URL="https://swift.rc.nectar.org.au:8888/v1/AUTH_7ea859948c3a451c9baced6fee813ed1/cpipe-assets-2.3"
MANIFEST="manifest-2.3"

# Helper functions
function err() {
    prefix=`date "+%F %T"`
    echo
    echo "========================= ERROR =================================="
    echo
    echo "$prefix $1" | fmt -w 100
    echo
    echo "=================================================================="
    echo
    exit 1
}

function warn() {
    prefix=`date "+%F %T"`
    echo
    echo "================================================================"
    echo "$prefix WARNING: $1" | fmt -w 100
    echo "================================================================"
    echo
}


function msg() {
    prefix=`date "+%F %T"`
    echo
    echo "================================================================"
    echo "$prefix $1"
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
    CONFIG=`sed 's/\/\/.*$//' $BASE/pipeline/config.groovy` 
    eval "$CONFIG"
}

function set_config_variable() {
    NAME="$1"
    VALUE="$2"
    cp "$BASE/pipeline/config.groovy" "$BASE/pipeline/config.groovy.tmp"
    sed 's,'^[\s]*$NAME'=\("\?\).*$,'$NAME'=\1'$VALUE'\1,g' $BASE/pipeline/config.groovy.tmp > "$BASE/pipeline/config.groovy" || err "Failed to set configuration variable $NAME to value $VALUE"
    rm "$BASE/pipeline/config.groovy.tmp"
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
            cp -v pipeline/config.groovy.template pipeline/config.groovy
        fi

        BASE=`pwd`
        set_config_variable BASE "$BASE"
        load_config
}



load_config

msg "Check base location $BASE is correct ..."
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

compile "$BCFTOOLS/bcftools"

compile "$HTSLIB/tabix"

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

########## vep ##########
msg "Check VEP database downloaded for version $VEP_VERSION..."
if [ -e $VEP/../vep_cache/homo_sapiens/${VEP_VERSION}*/1 ] && [ -e $VEP/Bio ]; then
    msg "VEP installed..."
else
    echo "
    Cpipe uses the Variant Effect Predictor from Ensembl
    to perform annotation of variants.
    
    To work, VEP must first be installed, and then reference
    files must be downloaded. Would you like to launch the
    VEP installer script now?

    NOTE: VEP may take several hours to download and install 
          reference data.
    "

    prompt "Do you want to run the VEP installer now? (y/n)" "y"
    if [ "$REPLY" == "y" ];
    then
        cd $VEP; 
        perl INSTALL.pl --CACHEDIR ../vep_cache --AUTO acf --SPECIES homo_sapiens_vep,homo_sapiens_refseq,homo_sapiens_merged --ASSEMBLY GRCh37 || err "Failed to run VEP installer"
        perl convert_cache.pl -species homo_sapiens -version ${VEP_VERSION}_GRCh37 --dir ../vep_cache || err "Failed to run VEP tabix for homo_sapiens"
        perl convert_cache.pl -species homo_sapiens_refseq -version ${VEP_VERSION}_GRCh37 --dir ../vep_cache || err "Failed to run VEP tabix for homo_sapiens_refseq"
        perl convert_cache.pl -species homo_sapiens_merged -version ${VEP_VERSION}_GRCh37 --dir ../vep_cache || err "Failed to run VEP tabix for homo_sapiens_merged"
    else
        msg "WARNING: Cpipe will not operate correctly if VEP is not installed"
    fi
fi

########## condel plugin ##########
msg "Configuring Condel Plugin ..."
cp "$CONDEL/config/condel_SP.conf.template" "$CONDEL/config/condel_SP.conf"

if [ ! -e $TOOLS/vep_plugins/Condel.pm ]; then
  ln -s "$CONDEL/Condel.pm" "$TOOLS/vep_plugins"
else
  msg "condel symlink already configured"
fi

sed -i 's,do not use,'$CONDEL/config',' $CONDEL/config/condel_SP.conf || err "Unable to configure Condel plugin"

########## dbnsfp plugin ##########
msg "Configuring dbNSFP plugin"

if [ ! -e $TOOLS/vep_plugins/dbNSFP.pm ]; then
  ln -s "$DBNSFP/dbNSFP.pm" "$TOOLS/vep_plugins"
else
  msg "condel symlink already configured"
fi

# download dbnsfp dataset
if [ ! -e "$TOOLS/vep_plugins/dbNSFP/dbNSFPv3.0b2a.zip" ]; then
  pushd "$TOOLS/vep_plugins/dbNSFP"
  DBNSFP_URL="ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.0b2a.zip"
  msg "downloading dbnsfp..."
  wget $DBNSFP_URL || err "Failed to download $DBNSFP_URL"
  msg "processing dbnsfp..."
  unzip dbNSFPv3.0b2a.zip
  cat dbNSFP*chr* | "$HTSLIB/bgzip" -c > dbNSFP.gz
  "$HTSLIB/tabix" -s 1 -b 2 -e 2 dbNSFP.gz
  msg "processing dbnsfp: done"
  popd
else
  msg "dbnsfp dataset already downloaded"
fi

##########
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
      if [ -x "$(command -v mysql)" ];
      then
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
            "select chrom, size from hg19.chromInfo"  > $HG19_CHROM_INFO || err "Failed to download genome file"
      else
        warn "No mysql client found. Cpipe will not operate correctly unless the $HG19_CHROM_INFO file is created"
        prompt "Press enter to continue" " "
      fi
    else
        warn "Cpipe will not operate correctly unless the $HG19_CHROM_INFO file is created"
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
    warn "The limit on open files is set to $MAX_OPEN_FILES. Cpipe may require more open files than this. Consider adding 'ulimit -S -n 2048' to your .bashrc file."
fi

msg "Check versions up to date ..."
if [ `basename $GROOVY_NGS` == "1.0" ];
then
    prompt "Your configuration is set to use an outdated version of groovy-ngs-utils. Update now (y/n)?" "y"
    if [ "$REPLY" == "y" ];
    then
        set_config_variable GROOVY_NGS "$TOOLS/groovy-ngs-utils/1.0.1"
    fi
fi

msg "Checking for custom assets..."
# download the manifest
pushd $REFBASE
[ -f $MANIFEST ] && rm $MANIFEST
wget $CUSTOM_ASSETS_URL/$MANIFEST

while read line; do
  if [[ "$line" =~ ^#.* ]]; then
    : #echo "comment"
  else
    filename=${line%,*}
    md5=${line##*,}
    echo "checking $filename..."
    if [ -f $filename ]; then
      # check md5
      existing=`md5sum $filename | awk '{ print $1; }'`
      if [ $existing == $md5 ]; then
        echo "$filename is up to date"
      else
        echo "updating $filename"
        rm $filename
        wget "$CUSTOM_ASSETS_URL/$filename"
      fi
    else # download file
      echo "downloading $filename"
      wget "$CUSTOM_ASSETS_URL/$filename"
    fi
  fi
done <"$MANIFEST"
popd 

msg "Success: all the dependencies I know how to check are OK"

