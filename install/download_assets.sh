#!/usr/bin/env bash

### Version constants ###
BWA_VERSION="0.7.13"
HTSLIB_VERSION="1.3" # Samtools and Bcftools also use this
BEDTOOLS_VERSION="2.25.0"
GATK_VERSION="3.6"
VEP_VERSION="84"
PYTHON_VERSION="2.7.12"
PERL_VERSION="5.24.0"
R_VERSION="3.3.1"
GROOVY_VERSION="2.4.7"
CPSUITE_VERSION="1.2.7"
FASTQC_VERSION="0.11.5"

ROOT=$( cd "$(dirname "${BASH_SOURCE}")"/.. ; pwd -P ) #The cpipe root directory
TOOLS_ROOT=$ROOT/tools
DATA_ROOT=$ROOT/data
PYTHON_ROOT=$TOOLS_ROOT/python
PERL_ROOT=$TOOLS_ROOT/perl
R_ROOT=$TOOLS_ROOT/r
JAVA_LIBS_ROOT=$TOOLS_ROOT/java_libs
GROOVY_ROOT=$TOOLS_ROOT/groovy

### Utility Functions ###
function download_gz {
    # $1 is the URL
    # $2 is the directory
    # e.g. download_gz http://cran.csiro.au/src/base/R-3/R-3.3.1.tar.gz /mnt/cpipe/tools/r

    FILE_NAME=`basename $1`\
    && mkdir -p $2\
    && curl -s $1 | tar -xz --strip-components=1 -C $2

    # && echo -n "Downloading $FILE_NAME into $2..."\
    check_success
}

function download_zip {
    # $1 is the URL
    # $2 is the directory to extract into
    # e.g. download_zip https://dl.bintray.com/groovy/maven/apache-groovy-binary-2.4.7.zip /mnt/cpipe/tools/groovy
    #      && echo -n "Downloading $FILE_NAME into $2..."\
      FILE_NAME=`basename $1` `#FILE_NAME is just the zip file without a URL e.g. apache-groovy-binary-2.4.7.zip`\
      && mkdir -p $2 `# Make the target directory`\
      && ZIP_FILE=$2/$FILE_NAME `# ZIP_FILE is the full path to the downloaded zip file`\
      && wget $1 -P $2 -q `#Perform the download`\
      && unzip -d $2 -qq $ZIP_FILE `#Unzip the file`\
      && rm $ZIP_FILE `#Delete the zip file`\
      && ZIP_ROOT=$2/`ls $2` `#Work out the root directory inside the zip file`\
      && mv $ZIP_ROOT/* $2 `#Move the contents out of this subdirectory because it doesn\'t have a consistent name`\
      && rm -r $ZIP_ROOT

      check_success
}

function check_success {
    if [ $? -ne 0 ] ; then
        echo failed.
    else
        echo success!
    fi
}

function existsExactlyOne {
# Fails unless there is exactly one argument which is a filename to an existing file
    [[ $# -eq 1 && -f $1 ]];
}

### Preliminary checks ###
## Check Java ##
if [ -f $JAVA_HOME/bin/java ] ; then
    JAVA_VER=`$JAVA_HOME/bin/java -version 2>&1 | sed -nr 's/.*version "([0-9]+\.[0-9]+).*/\1/p'`
    VALID_JAVA=`awk -v ver=$JAVA_VER 'BEGIN{ print (ver < 1.8) ? "FAILED" : "PASSED" }'`
else
    JAVA_VER=NONE
    VALID_JAVA=FAILED
fi

if [ $VALID_JAVA != 'PASSED' ]; then
    echo 'Your java version according to JAVA_HOME is' $JAVA_VER'. Please install Java 1.8 or greater to compile GATK'
    exit 1
fi

function pushd {
    command pushd "$@" > /dev/null
}

function popd {
    command popd "$@" > /dev/null
}

### Start of script ###

## General Dependencies ##
sudo cpan App::cpanminus > /dev/null
export JAVA_HOME=/usr

## Move Paths ##
mkdir -p $DATA_ROOT $TOOLS_ROOT $JAVA_LIBS_ROOT
cd $TOOLS_ROOT

##Language installations##
#Python
if [[ ! -e $PYTHON_ROOT ]]; then
    download_gz https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz $PYTHON_ROOT
fi

#Perl
if [[ ! -e $PERL_ROOT ]]; then
    download_gz http://www.cpan.org/src/5.0/perl-$PERL_VERSION.tar.gz $PERL_ROOT\
      && mv $PERL_ROOT/configure.gnu $PERL_ROOT/configure.sh
fi

#R
if [[ ! -e $R_ROOT ]]; then
    download_gz http://cran.csiro.au/src/base/R-3/R-$R_VERSION.tar.gz $R_ROOT
fi

#Groovy
if [[ ! -e $GROOVY_ROOT ]]; then
    download_zip https://dl.bintray.com/groovy/maven/apache-groovy-binary-$GROOVY_VERSION.zip $TOOLS_ROOT/groovy\
    && mv $TOOLS_ROOT/groovy-$GROOVY_VERSION $TOOLS_ROOT/groovy
fi

#BWA
if [[ ! -e $TOOLS_ROOT/bwa ]]; then
    echo -n 'Downloading BWA...'\
    && download_gz https://codeload.github.com/lh3/bwa/tar.gz/v$BWA_VERSION $TOOLS_ROOT/bwa
    check_success
fi

#Htslib (requirement for samtools and bcftool)
if [[ ! -e $TOOLS_ROOT/htslib ]]; then
    echo -n 'Downloading Htslib...'\
     && download_gz https://codeload.github.com/samtools/htslib/tar.gz/HTSLIB_VERSION $TOOLS_ROOT/htslib
    check_success
fi

#Samtools
if [[ ! -e $TOOLS_ROOT/samtools ]]; then
    echo -n 'Downloading Samtools...'\
    && download_gz https://codeload.github.com/samtools/samtools/tar.gz/HTSLIB_VERSION $TOOLS_ROOT/samtools
    check_success
fi

#Bcftools
if [[ ! -e $TOOLS_ROOT/bwa ]]; then
    echo -n 'Downloading Bcftools...'\
    && download_gz https://codeload.github.com/samtools/bcftools/tar.gz/HTSLIB_VERSION $TOOLS_ROOT/bcftools
    check_success
fi

#Bedtools
if [[ ! -e $TOOLS_ROOT/bedtools ]]; then
    echo -n 'Downloading Bedtools...'\
    && download_gz https://codeload.github.com/arq5x/bedtools2/tar.gz/v$BEDTOOLS_VERSION $TOOLS_ROOT/bedtools
    check_success
fi

#VEP, including assets. Involves downloading ensembl-tools and deleting everything that isn't the VEP script
if [[ ! -e $TOOLS_ROOT/vep ]]; then
    echo -n 'Downloading VEP...'\
        && wget https://github.com/Ensembl/ensembl-tools/archive/release/$VEP_VERSION.zip -q -O vep.zip\
        && unzip vep.zip > /dev/null\
        && mv ensembl-tools-release-$VEP_VERSION/scripts/variant_effect_predictor vep\
        && rm -rf ensembl-tools-release-$VEP_VERSION vep.zip
    check_success
fi

# Fastqc
if [[ ! -e $TOOLS_ROOT/fastqc ]]; then
    echo -n 'Downloading fastqc...'\
    && download_zip "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" $TOOLS_ROOT/fastqc
    check_success
fi

# Bpipe
if [[ ! -e $TOOLS_ROOT/bpipe ]]; then
    git clone https://github.com/ssadedin/bpipe
    pushd bpipe
        ./gradlew dist
    popd
fi

# Setup Perl variables
PERL5LIB=$TOOLS_ROOT:$PERL5LIB
PATH=$TOOLSROOT/htslib:$PATH

echo -n 'Installing VEP perl dependencies...'\
&& pushd $TOOLS_ROOT/vep\
    && { [[ -e $TOOLS_ROOT/vep/cpanfile ]] || mv $ROOT/cpanfile $TOOLS_ROOT/vep ;}\
    && sudo cpanm --installdeps . > /dev/null\
&& popd
check_success

## Data Files ##

if [[ ! -e $DATA_ROOT/vep_cache ]]; then
    echo -n 'Installing VEP databases...'\
    && VEP_CACHE=$DATA_ROOT/vep_cache\
    && mkdir $VEP_CACHE\
    && perl $TOOLS_ROOT/vep/INSTALL.pl --CACHEDIR $VEP_CACHE --AUTO acf --SPECIES homo_sapiens_vep,homo_sapiens_refseq,homo_sapiens_merged --ASSEMBLY GRCh37 > /dev/null
    check_success
fi

if [[ ! -e $DATA_ROOT/gatk ]]; then
    echo -n 'Downloading GATK bundle...'
    mkdir $DATA_ROOT/gatk
    GATK_BUNDLE_ROOT=ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
    GATK_BUNDLE_FILES="dbsnp_138.hg19.vcf.gz\
    dbsnp_138.hg19.vcf.idx.gz\
    Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz\
    Mills_and_1000G_gold_standard.indels.hg19.vcf.idx.gz\
    ucsc.hg19.dict.gz\
    ucsc.hg19.fasta.gz\
    ucsc.hg19.fasta.fai.gz"

    for f in $GATK_BUNDLE_FILES ;  do
         URL="$GATK_BUNDLE_ROOT$f";
         BASE=`basename $f .gz`;
         echo "Downloading $BASE...";
         curl --user gsapubftp-anonymous:cpipe.user@cpipeline.org $URL | gunzip > $DATA_ROOT/gatk/$BASE;
         check_success;
    done
fi

if [[ ! -f $DATA_ROOT/hg19.genome ]]; then
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" > $DATA_ROOT/hg19.genome
fi

## Index Reference File ##
$TOOLS_ROOT/bwa/bwa index -a bwtsw $DATA_ROOT/gatk/ucsc.hg19.fasta

## Jar Dependencies ##
pushd $JAVA_LIBS_ROOT
    #GATK, also pre-compile the .jar file
     if ! existsExactlyOne $JAVA_LIBS_ROOT/GenomeAnalysisTK.jar ; then
        echo -n 'Downloading GATK...'\
        && download_gz https://codeload.github.com/broadgsa/gatk-protected/tar.gz/$GATK_VERSION $JAVA_LIBS_ROOT/gatk
        check_success

        echo -n 'Compiling GATK...'
        pushd gatk\
            && mvn --quiet verify -P\!queue > /dev/null\
            && mv target/executable/GenomeAnalysisTK.jar $JAVA_LIBS_ROOT\
            && shopt -s extglob\
        && popd\
        && rm -rf gatk
        check_success
    fi

    if ! existsExactlyOne $JAVA_LIBS_ROOT/JUnitXmlFormatter*.jar ; then
        echo "Compiling JUnitXmlFormatter dependencies"\
        && git clone https://github.com/barrypitman/JUnitXmlFormatter\
        && pushd JUnitXmlFormatter\
            && mvn install\
            && mv target/JUnitXmlFormatter* $JAVA_LIBS_ROOT\
        && popd\
        && rm -rf JUnitXmlFormatter
        check_success
    fi

    # Groovy ngs utils
    if ! existsExactlyOne $JAVA_LIBS_ROOT/groovy-ngs-utils.jar ; then
        echo "Compiling groovy-ngs-utils dependencies"\
        && git clone https://github.com/ssadedin/groovy-ngs-utils -b upgrade-biojava --depth=1 --quiet\
        && pushd groovy-ngs-utils\
        && ./gradlew jar --quiet\
        && popd\
        && mv $JAVA_LIBS_ROOT/groovy-ngs-utils/build/libs/groovy-ngs-utils.jar $JAVA_LIBS_ROOT\
        && rm -rf groovy-ngs-utils
    fi


    if ! existsExactlyOne $JAVA_LIBS_ROOT/takari-cpsuite* ; then
        echo "Downloading cpsuite"
        mvn --quiet dependency:copy \
            -Dartifact=io.takari.junit:takari-cpsuite:$CPSUITE_VERSION\
            -DoutputDirectory=$JAVA_LIBS_ROOT\
            -DstripVersion=true
        check_success
    fi

popd

# Note: this has been commented out as tabix format does not seem to store frequency information
#echo -n 'Converting cache to tabix format...'
#perl convert_cache.pl -species homo_sapiens -version ${VEP_VERSION}_GRCh37 --dir $VEP_CACHE || echo "Failed to run VEP tabix for homo_sapiens"\
#    && perl convert_cache.pl -species homo_sapiens_refseq -version ${VEP_VERSION}_GRCh37 --dir $VEP_CACHE || echo "Failed to run VEP tabix for homo_sapiens_refseq"\
#    && perl convert_cache.pl -species homo_sapiens_merged -version ${VEP_VERSION}_GRCh37 --dir $VEP_CACHE || echo "Failed to run VEP tabix for homo_sapiens_merged"\
#check_success
cd ..
