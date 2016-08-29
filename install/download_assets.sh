#!/usr/bin/env bash

### Imports ###
source `dirname ${BASH_SOURCE[0]}`/common/compile.sh

### Logging ####
LOG_FILE=$(readlink -f $(dirname ${BASH_SOURCE[0]})/download.log)
echo "## Running downloader on `date` ##" > $LOG_FILE
echo "## Logging to $LOG_FILE ##"

### Version constants ###
BWA_VERSION="0.7.13"
HTSLIB_VERSION="1.3" # Samtools and Bcftools also use this
BEDTOOLS_VERSION="2.25.0"
GATK_VERSION="3.6"
VEP_VERSION="85"
PYTHON_VERSION="2.7.12"
PERL_VERSION="5.24.0"
R_VERSION="3.3.1"
GROOVY_VERSION="2.4.7"
CPSUITE_VERSION="1.2.7"
FASTQC_VERSION="0.11.5"
PICARD_VERSION="2.6.0"

ROOT=$(readlink -f $(dirname ${BASH_SOURCE[0]})/..) #The cpipe root directory
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
    {
        FILE_NAME=`basename $1`\
        && mkdir -p $2\
        && curl $1 2>> $LOG_FILE| tar -xz --strip-components=1 -C $2
    } >> $LOG_FILE
}

function download_zip {
    # $1 is the URL
    # $2 is the directory to extract into
    # e.g. download_zip https://dl.bintray.com/groovy/maven/apache-groovy-binary-2.4.7.zip /mnt/cpipe/tools/groovy
    #      && echo -n "Downloading $FILE_NAME into $2..."\
   {
        FILE_NAME=`basename $1` `#FILE_NAME is just the zip file without a URL e.g. apache-groovy-binary-2.4.7.zip`\
        && mkdir -p $2 `# Make the target directory`\
        && ZIP_FILE=$2/$FILE_NAME `# ZIP_FILE is the full path to the downloaded zip file`\
        && wget $1 -P $2 `#Perform the download`\
        && unzip -d $2 $ZIP_FILE `#Unzip the file`\
        && rm $ZIP_FILE `#Delete the zip file`\
        && ZIP_ROOT=$2/`ls $2` `#Work out the root directory inside the zip file`\
        && mv $ZIP_ROOT/* $2 `#Move the contents out of this subdirectory because it doesn\'t have a consistent name`\
        && rm -r $ZIP_ROOT
    } >> $LOG_FILE
}

function check_success {
    if [ $? -ne 0 ] ; then
        echo failed.
        echo "Look at $LOG_FILE for details"
        exit 1
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

function command_exists {
    type $1 > /dev/null 2>&1
}

### Start of script ###

{
    ## Move Paths ##
    mkdir -p $DATA_ROOT $TOOLS_ROOT $JAVA_LIBS_ROOT
    pushd $TOOLS_ROOT

        ## General Dependencies ##
        #echo -n 'Installing cpanm...'
        #if command_exists cpanm; then
        #    echo 'already satisfied.'
        #else
        #    sudo cpan App::cpanminus >> $LOG_FILE
        #    check_success
        #fi

        ##Language installations##
        #Python
        echo -n 'Downloading python...'
        if [[ ! -e $PYTHON_ROOT ]]; then
            download_gz https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz $PYTHON_ROOT
            check_success
        else
            echo 'already satisfied'
        fi


        #Perl
        echo -n 'Downloading perl...'
        if [[ ! -e $PERL_ROOT ]]; then
            download_gz http://www.cpan.org/src/5.0/perl-$PERL_VERSION.tar.gz $PERL_ROOT\
            && mv $PERL_ROOT/configure.gnu $PERL_ROOT/configure.sh
            check_success
        else
            echo 'already satisfied'
        fi

        #R
        echo -n 'Downloading R...'
        if [[ ! -e $R_ROOT ]]; then
            download_gz http://cran.csiro.au/src/base/R-3/R-$R_VERSION.tar.gz $R_ROOT
            check_success
        else
            echo 'already satisfied'
        fi

        #Groovy
        echo -n 'Downloading groovy...'
        if [[ ! -e $GROOVY_ROOT ]]; then
            download_zip https://dl.bintray.com/groovy/maven/apache-groovy-binary-$GROOVY_VERSION.zip $TOOLS_ROOT/groovy
            check_success
        else
            echo 'already satisfied'
        fi

        #BWA. Also compile in order to perform indexing
        echo -n 'Downloading BWA...'
        if [[ ! -e $TOOLS_ROOT/bwa ]]; then
            download_gz https://codeload.github.com/lh3/bwa/tar.gz/v$BWA_VERSION $TOOLS_ROOT/bwa
            check_success

            echo -n 'Compiling BWA...'
            compile $TOOLS_ROOT/bwa >> $LOG_FILE
            check_success
        else
            echo 'already satisfied'
        fi

        #Htslib (requirement for samtools and bcftool)
        echo -n 'Downloading Htslib...'
        if [[ ! -e $TOOLS_ROOT/htslib ]]; then
            download_gz https://codeload.github.com/samtools/htslib/tar.gz/$HTSLIB_VERSION $TOOLS_ROOT/htslib
            check_success
        else
            echo 'already satisfied'
        fi

        #Samtools. Also compile in order to perform indexing
        echo -n 'Downloading Samtools...'
        if [[ ! -e $TOOLS_ROOT/samtools ]]; then
            download_gz https://codeload.github.com/samtools/samtools/tar.gz/$HTSLIB_VERSION $TOOLS_ROOT/samtools
            check_success

             echo -n 'Compiling Samtools...'
             compile $TOOLS_ROOT/samtools >> $LOG_FILE
             check_success
        else
            echo 'already satisfied'
        fi

        #Bcftools
        echo -n 'Downloading Bcftools...'
        if [[ ! -e $TOOLS_ROOT/bcftools ]]; then
            download_gz https://codeload.github.com/samtools/bcftools/tar.gz/$HTSLIB_VERSION $TOOLS_ROOT/bcftools
            check_success
        else
            echo 'already satisfied'
        fi

        #Bedtools
        echo -n 'Downloading Bedtools...'
        if [[ ! -e $TOOLS_ROOT/bedtools ]]; then
            download_gz https://codeload.github.com/arq5x/bedtools2/tar.gz/v$BEDTOOLS_VERSION $TOOLS_ROOT/bedtools
            check_success
        else
            echo 'already satisfied'
        fi

        #VEP, including assets. Involves downloading ensembl-tools and deleting everything that isn't the VEP script
        echo -n 'Downloading VEP...'
        if [[ ! -e $TOOLS_ROOT/vep ]]; then
            wget https://github.com/Ensembl/ensembl-tools/archive/release/$VEP_VERSION.zip -O vep.zip >> $LOG_FILE\
            && unzip vep.zip >> $LOG_FILE\
            && mv ensembl-tools-release-$VEP_VERSION/scripts/variant_effect_predictor vep\
            && rm -rf ensembl-tools-release-$VEP_VERSION vep.zip
            check_success
        else
            echo 'already satisfied'
        fi

        # Fastqc
        echo -n 'Downloading fastqc...'
        if [[ ! -e $TOOLS_ROOT/fastqc ]]; then
            download_zip "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" $TOOLS_ROOT/fastqc
            check_success
        else
            echo 'already satisfied'
        fi

        # Bpipe
        echo -n 'Downloading bpipe...'
        if [[ ! -e $TOOLS_ROOT/bpipe ]]; then
            git clone https://github.com/ssadedin/bpipe\
            && pushd bpipe\
                && ./gradlew dist >> $LOG_FILE\
            && popd
            check_success
        else
            echo 'already satisfied'
        fi

        #GATK, also pre-compile the .jar file
        echo -n 'Downloading GATK...'
        if [[ ! -e $TOOLS_ROOT/gatk ]]; then
            mkdir -p $TOOLS_ROOT/gatk\
            && download_gz https://codeload.github.com/broadgsa/gatk-protected/tar.gz/$GATK_VERSION $TOOLS_ROOT/gatk
            check_success

            echo -n 'Compiling GATK...'
            pushd $TOOLS_ROOT/gatk\
                && mvn verify -P\!queue >> $LOG_FILE\
                && GATK_JAR=`readlink -f target/GenomeAnalysisTK.jar`\
                && unlink target/GenomeAnalysisTK.jar\
                && mv $GATK_JAR ./GenomeAnalysisTK.jar\
                && bash -O extglob -c 'rm -rf !(GenomeAnalysisTK.jar)'\
            && popd
            check_success
        else
            echo 'already satisfied'
        fi

        echo -n 'Downloading picard...'
        if [[ ! -e $TOOLS_ROOT/picard ]]; then
            mkdir -p $TOOLS_ROOT/picard\
            && wget -P $TOOLS_ROOT/picard https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar >> $LOG_FILE
            check_success
        else
            echo 'already satisfied'
        fi

        # Setup Perl variables
        export PERL5LIB=$TOOLS_ROOT/perl_lib/lib/perl5:$PERL5LIB
        export PATH=$TOOLSROOT/htslib:$PATH

        echo -n 'Installing general perl dependencies...'
         if [[ ! -e $TOOLS_ROOT/perl_lib ]] ; then
            mkdir -p $TOOLS_ROOT/perl_lib\
            && pushd $ROOT/install\
                && cpanm --installdeps --local-lib-contained $TOOLS_ROOT/perl_lib . >> $LOG_FILE\
            && popd
            check_success
        else
             echo 'already done'
        fi

        echo -n 'Installing VEP dependencies ...'
        if [[ ! -e $TOOLS_ROOT/perl_lib/lib/perl5/Bio/EnsEMBL ]]; then
            pushd $TOOLS_ROOT/vep\
                && yes | perl $TOOLS_ROOT/vep/INSTALL.pl --NO_HTSLIB --AUTO a --DESTDIR $TOOLS_ROOT/perl_lib/lib/perl5  >> $LOG_FILE\
            && popd
            check_success
        else
             echo 'already done'
        fi

        echo -n 'Installing VEP plugins...'
        if [[ ! -f $TOOLS_ROOT/perl_lib/lib/perl5/Condel.pm ]] ; then
            git clone https://github.com/Ensembl/VEP_plugins\
            && mv VEP_plugins/*.pm $TOOLS_ROOT/perl_lib/lib/perl5\
            && rm -rf VEP_plugins
            check_success
        else
             echo 'already done'
        fi

        ## Data Files ##
        echo -n 'Installing VEP cache and reference file...'
        if [[ ! -e $DATA_ROOT/vep_cache ]]; then
            # Note that if you include more than 1 species then the assembly fasta file will only be installed into the last
            VEP_CACHE=$DATA_ROOT/vep_cache\
            && mkdir $VEP_CACHE\
            && perl $TOOLS_ROOT/vep/INSTALL.pl --NO_HTSLIB --CACHEDIR $VEP_CACHE --AUTO cf --SPECIES homo_sapiens_refseq --ASSEMBLY GRCh37  >> $LOG_FILE
            check_success
        else
             echo 'already done'
        fi

        echo 'Downloading GATK bundle...'
        if [[ ! -e $DATA_ROOT/gatk ]]; then
            mkdir $DATA_ROOT/gatk
        fi

        GATK_BUNDLE_ROOT=ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
        GATK_BUNDLE_FILES="dbsnp_138.hg19.vcf.gz\
        dbsnp_138.hg19.vcf.idx.gz\
        Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz\
        Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz\
        ucsc.hg19.dict.gz\
        ucsc.hg19.fasta.gz\
        ucsc.hg19.fasta.fai.gz"

        for f in $GATK_BUNDLE_FILES ;  do
             URL="$GATK_BUNDLE_ROOT$f"
             BASE=`basename $f .gz`
             echo -e -n "\tDownloading $BASE..."
             if [[ ! -e $DATA_ROOT/gatk/$BASE ]]; then
                 curl --user gsapubftp-anonymous:cpipe.user@cpipeline.org $URL | gunzip > $DATA_ROOT/gatk/$BASE
                 check_success
             else
                echo "already exists"
             fi
        done

        echo -n 'Downloading chromosome sizes...'
        if [[ ! -f $DATA_ROOT/chromosomes/hg19.genome ]]; then
            mkdir -p $DATA_ROOT/chromosomes\
            && mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" > $DATA_ROOT/chromosomes/hg19.genome
            check_success
        else
            echo "already exists"
        fi

        ## Index Reference File ##
        echo -n 'Indexing reference file using bwa...'
        if [[ ! -f $DATA_ROOT/gatk/ucsc.hg19.fasta.bwt ]] ; then
            $TOOLS_ROOT/bwa/bwa index -a bwtsw $DATA_ROOT/gatk/ucsc.hg19.fasta
            check_success
        else
            echo "already done"
        fi

        ## Index Reference File ##
        echo -n 'Indexing reference file using samtools...'
        if [[ ! -f $DATA_ROOT/gatk/ucsc.hg19.fasta.fai ]] ; then
            $TOOLS_ROOT/samtools/samtools faidx $DATA_ROOT/gatk/ucsc.hg19.fasta
            check_success
        else
            echo "already done"
        fi

        ## Jar Dependencies ##
        pushd $JAVA_LIBS_ROOT
            echo -n "Downloading and compiling JUnitXmlFormatter..."
            if ! existsExactlyOne $JAVA_LIBS_ROOT/JUnitXmlFormatter*.jar ; then
                git clone https://github.com/barrypitman/JUnitXmlFormatter\
                && pushd JUnitXmlFormatter\
                    && mvn install >> $LOG_FILE\
                    && mv target/JUnitXmlFormatter* $JAVA_LIBS_ROOT\
                && popd\
                && rm -rf JUnitXmlFormatter
                check_success
            else
                echo "already done"
            fi

            # Groovy ngs utils
            echo -n "Downloading and compiling groovy-ngs-utils..."
            if ! existsExactlyOne $JAVA_LIBS_ROOT/groovy-ngs-utils.jar ; then
                git clone https://github.com/ssadedin/groovy-ngs-utils -b upgrade-biojava --depth=1 --quiet\
                && pushd groovy-ngs-utils\
                && ./gradlew jar >> $LOG_FILE\
                && popd\
                && mv $JAVA_LIBS_ROOT/groovy-ngs-utils/build/libs/groovy-ngs-utils.jar $JAVA_LIBS_ROOT\
                && rm -rf groovy-ngs-utils
                check_success
            else
                echo "already done"
            fi

            echo -n "Downloading and compiling takari-cpsuite..."
            if ! existsExactlyOne $JAVA_LIBS_ROOT/takari-cpsuite* ; then
                echo "Downloading cpsuite"
                mvn dependency:copy \
                    -Dartifact=io.takari.junit:takari-cpsuite:$CPSUITE_VERSION\
                    -DoutputDirectory=$JAVA_LIBS_ROOT\
                    -DstripVersion=true >> $LOG_FILE
                check_success
            else
                echo "already done"
            fi

        popd
    popd
} 2>> $LOG_FILE # Redirect all errors to the log file

echo "Done!"
