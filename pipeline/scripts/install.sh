#!/bin/bash
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

# 
# Load config from config.groovy
# If an argument is passed then the values are filtered through grep to only
# include the given values.
#
# Note that here we are taking advantage of the fact that Groovy and Bash share
# common syntax for defining variables.
#
function load_config() {
    if [ -z "$BASE" ];
    then
        BASE="."
    fi
    
    if [ ! -z "$1" ];
    then
        CONFIG=`sed 's/\/\/.*$//' $BASE/pipeline/config.groovy | grep "$1"` 
    else
        CONFIG=`sed 's/\/\/.*$//' $BASE/pipeline/config.groovy` 
    fi
    
    eval "$CONFIG"
}

function set_config_variable() {
    NAME="$1"
    VALUE="$2"
    cp "$BASE/pipeline/config.groovy" "$BASE/pipeline/config.groovy.tmp"
   
    $GROOVY -D name="$NAME" -D value="$VALUE" \
      -pne 'line.startsWith(System.properties.name+"=")?line.replaceAll("=.*",/="/+java.util.regex.Matcher.quoteReplacement(System.properties.value)+/"/): line' \
      "$BASE/pipeline/config.groovy.tmp" > \
      "$BASE/pipeline/config.groovy" \
        || err "Failed to set configuration variable $NAME to value $VALUE"
        
    rm "$BASE/pipeline/config.groovy.tmp"
    load_config
}

function gatk_prompt() {
        echo "
 The recommended version of GATK for Cpipe is 3.5. However license
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
            set_config_variable GATK "$TOOLS/gatk/2.3.9"
            set_config_variable GATK_LEGACY "true"
        else
            msg "WARNING: your installation will not work unless you set GATK manually youself in pipeline/config.groovy"
        fi
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
        TOOLS="$BASE/tools"
        load_config GROOVY
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
    if [ -e $TOOLS/gatk/3.5 ];
    then    
       prompt "Found GATK 3.5. Continue with this version? (y/n)" "y"
       if [ "$REPLY" == "y" ];
       then
           set_config_variable GATK '$TOOLS/gatk/3.5'
       else
           gatk_prompt
       fi    
    else
        err "Please install GATK then re-run this script"
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
        cd $VEP
        export PERL5LIB="$PERL5LIB:$TOOLS/perl5:$TOOLS/perl5/lib/perl5"
        # convert ../vep_cache to absolute path
        VEP_CACHE=`echo "$VEP" | sed 's/\/[^\/]*$/\/vep_cache/'`
        msg "INFO: VEP is installing homo_sapiens_vep"
        perl INSTALL.pl --CACHEDIR $VEP_CACHE --AUTO acf --SPECIES homo_sapiens_vep --ASSEMBLY GRCh37 || err "Failed to run VEP installer"
        msg "INFO: VEP is installing homo_sapiens_refseq"
        perl INSTALL.pl --CACHEDIR $VEP_CACHE --AUTO acf --SPECIES homo_sapiens_refseq --ASSEMBLY GRCh37 || err "Failed to run VEP installer"
        msg "INFO: VEP is installing homo_sapiens_merged"
        perl INSTALL.pl --CACHEDIR $VEP_CACHE --AUTO acf --SPECIES homo_sapiens_merged --ASSEMBLY GRCh37 || err "Failed to run VEP installer"
        # we don't run convert_cache as it (currently) messes up the frequency data (gmaf, etc)
        cd -
    else
        msg "WARNING: Cpipe will not operate correctly if VEP is not installed"
    fi
    msg "INFO: VEP has finished"
fi

########## condel plugin ##########
msg "Configuring Condel plugin..."
cp "$CONDEL/config/condel_SP.conf.template" "$CONDEL/config/condel_SP.conf"

if [ ! -e $TOOLS/vep_plugins/Condel.pm ]; then
  ln -s "$CONDEL/Condel.pm" "$TOOLS/vep_plugins"
else
  msg "condel symlink already configured"
fi

# Used to use sed to set this as below, but not all versions of sed
# support alternate pattern separators (s,foo,bar,g) which makes it tricky to use
# for file paths - instead use some inline groovy to do it
unset GROOVY_HOME
./tools/groovy/2.3.4/bin/groovy -e 'new File(args[0]).text = new File(args[0]).text.replaceAll("do not use", args[1])' \
           $CONDEL/config/condel_SP.conf $CONDEL\/config \
           || err "Unable to configure Condel plugin"

########## dbnsfp plugin ##########
msg "Configuring dbNSFP plugin"

if [ ! -e $TOOLS/vep_plugins/dbNSFP.pm ]; then
  ln -s "$DBNSFP/dbNSFP.pm" "$TOOLS/vep_plugins"
else
  msg "condel symlink already configured"
fi

# download dbnsfp dataset
#if [ ! -e "$TOOLS/vep_plugins/dbNSFP/dbNSFPv3.0b2a.zip" ]; then
if [ ! -e "$TOOLS/vep_plugins/dbNSFP/dbNSFPv2.9.1.zip" ]; then
  pushd "$TOOLS/vep_plugins/dbNSFP"
  #DBNSFP_URL="ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.0b2a.zip"
  DBNSFP_URL="ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv2.9.1.zip"
  msg "downloading dbnsfp..."
  wget $DBNSFP_URL || err "Failed to download $DBNSFP_URL"
  msg "downloading dbnsfp: done"
  popd
else
  msg "dbnsfp dataset already downloaded"
fi

# process dbnsfp dataset
if [ ! -e "$TOOLS/vep_plugins/dbNSFP/dbNSFP.gz" ]; then
  pushd "$TOOLS/vep_plugins/dbNSFP"
  msg "processing dbnsfp..."
  unzip dbNSFPv2.9.1.zip
  cat dbNSFP*chr* | "$HTSLIB/bgzip" -c > dbNSFP.gz
  "$HTSLIB/tabix" -s 1 -b 2 -e 2 dbNSFP.gz
  msg "processing dbnsfp: done"
  popd
else
  msg "dbnsfp dataset already downloaded"
fi

# grantham plugin
msg "Configuring Grantham plugin..."
if [ ! -e "$TOOLS/vep_plugins/Grantham.pm" ]; then
  ln -s "$TOOLS/vep_plugins/grantham/20160614/Grantham.pm" "$TOOLS/vep_plugins"
else
  msg "grantham symlink already configured"
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

# Since Mac's don't ship with md5sum 
: ${MD5SUM:=md5sum}
type $MD5SUM > /dev/null 2>&1 || {
  type md5sum-lite > /dev/null 2>&1 || err "Cannot find a suitable md5sum utility. Please set MD5SUM environment variable"
  MD5SUM=md5sum-lite
  echo "Using MD5SUM=$MD5SUM"
}

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
      existing=`$MD5SUM $filename | awk '{ print $1; }'`
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

