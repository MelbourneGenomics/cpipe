# Version constants
BWA_VERSION="0.7.13"
HTSLIB_VERSION="1.3" # Samtools and Bcftools also use this
BEDTOOLS_VERSION="2.25.0"
GATK_VERSION="3.6"
VEP_VERSION="84"

ROOT=$(dirname $(pwd)) #The cpipe root directory
TOOLS_ROOT=$ROOT/tools
DATA_ROOT=$ROOT/data

# Utility functions
function download_github {
    # Downloads a gzip file of the source for a particular github release
    # This could be adapted to use git clone but for now it directly downloads the release
    # $1 should be author/project e.g. "broadgsa/gatk-protected"
    # $2 can be any string corresponding to a release, e.g. "v2.5"
    # Returns (echoes) the directory, so it can be used to set a config variable

    lib=`echo $1 | cut -f2 -d'/'`
    url=https://codeload.github.com/$1/tar.gz/$2
    dir=$lib-${2//'v'}

    { curl -sS $url | tar -xz; }\
    	&& mv $dir $lib >/dev/null

    echo $lib
}

# $1 is the URL, $2 is the name of the asset, $3 is the directory
# e.g. download_zip_asset http://apache.mirror.amaze.com.au/groovy/2.4.7/sources/apache-groovy-src-2.4.7.zip groovy
function download_zip_asset {
    ZIP_FILE="$ROOT/$3/$2.zip"\
    && wget $1 -O $ZIP_FILE\
    && unzip $ZIP_FILE -d "$ROOT/$3/$2"\
    && rm ZIP_FILE
}

function check_success {
    if [ $? -ne 0 ] ; then
        echo failed.
    else
        echo success!
    fi
}

# Check versions
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

##Start of script##
# Dependencies
sudo cpan App::cpanminus
export JAVA_HOME=/usr

#Paths
mkdir $DATA_ROOT $TOOLS_ROOT
cd $TOOLS_ROOT

#BWA#
echo -n 'Downloading BWA...'
download_github lh3/bwa v$BWA_VERSION > /dev/null
check_success

#Htslib (requirement for samtools and bcftool)#
echo -n 'Downloading Htslib...'
download_github samtools/htslib $HTSLIB_VERSION > /dev/null
check_success

#Samtools#
echo -n 'Downloading Samtools...'
download_github samtools/samtools $HTSLIB_VERSION > /dev/null
check_success

#Bcftools#
echo -n 'Downloading Bcftools...'
download_github samtools/bcftools $HTSLIB_VERSION > /dev/null
check_success

#Bedtools#
echo -n 'Downloading Bedtools...'
download_github arq5x/bedtools2 v$BEDTOOLS_VERSION >/dev/null
check_success

#GATK, also pre-compile the .jar file
echo -n 'Downloading GATK...'
GATK=`download_github broadgsa/gatk-protected $GATK_VERSION`\
    && mv $GATK gatk
check_success

echo -n 'Compiling GATK...'
pushd gatk\
    && mvn --quiet verify -P\!queue > /dev/null\
    && mv target/executable/GenomeAnalysisTK.jar .\
    && shopt -s extglob\
    && rm -rf !(GenomeAnalysisTK.jar)
check_success
popd

#VEP, including assets#
# Download and extract
echo -n 'Downloading VEP...'
wget https://github.com/Ensembl/ensembl-tools/archive/release/$VEP_VERSION.zip -q -O vep.zip\
    && unzip vep.zip > /dev/null\
    && mv ensembl-tools-release-$VEP_VERSION/scripts/variant_effect_predictor vep\
    && rm -rf ensembl-tools-release-$VEP_VERSION vep.zip
check_success

# Setup Perl variables
PERL5LIB=$TOOLS_ROOT:$PERL5LIB
PATH=$TOOLSROOT/htslib:$PATH

echo -n 'Installing VEP perl dependencies...'
cd $TOOLS_ROOT/vep\
    && mv $ROOT/cpanfile $TOOLS_ROOT/vep\
    && sudo cpanm --installdeps . > /dev/null
check_success

echo -n 'Installing VEP databases...'
VEP_CACHE=$DATA_ROOT/vep_cache
mkdir $VEP_CACHE/
    && perl $TOOLS_ROOT/vep/INSTALL.pl --CACHEDIR $VEP_CACHE --AUTO acf --SPECIES homo_sapiens_vep,homo_sapiens_refseq,homo_sapiens_merged --ASSEMBLY GRCh37 > /dev/null
check_success

cd ..

echo -n 'Downloading GATK bundle...'
mkdir $DATA_ROOT/gatk
GATK_BUNDLE_ROOT=ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
GATK_BUNDLE_FILES="dbsnp_138.hg19.vcf.gz Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz ucsc.hg19.dict.gz ucsc.hg19.fasta.gz ucsc.hg19.fasta.fai.gz"

for f in $GATK_BUNDLE_FILES ;  do
     URL="$GATK_BUNDLE_ROOT$f";
     BASE=`basename $f .gz`;
     echo "Downloading $BASE...";
     curl --user gsapubftp-anonymous:cpipe.user@cpipeline.org $URL | gunzip > $DATA_ROOT/gatk/$BASE;
     check_success;
 done

# Note: this has been commented out as tabix format does not seem to store frequency information
#echo -n 'Converting cache to tabix format...'
#perl convert_cache.pl -species homo_sapiens -version ${VEP_VERSION}_GRCh37 --dir $VEP_CACHE || echo "Failed to run VEP tabix for homo_sapiens"\
#    && perl convert_cache.pl -species homo_sapiens_refseq -version ${VEP_VERSION}_GRCh37 --dir $VEP_CACHE || echo "Failed to run VEP tabix for homo_sapiens_refseq"\
#    && perl convert_cache.pl -species homo_sapiens_merged -version ${VEP_VERSION}_GRCh37 --dir $VEP_CACHE || echo "Failed to run VEP tabix for homo_sapiens_merged"\
#check_success
cd ..
