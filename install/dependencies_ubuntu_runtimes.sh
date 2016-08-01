PYTHON_VERSION=2.7.12
PERL_VERSION=5.24.0
R_VERSION=3.3.1
GROOVY_VERSION=2.4.7
JAVA_VERSON=8u91-b14-0ubuntu4~16.04.1

ROOT=$(dirname $(pwd)) #The cpipe root directory
TOOLS_ROOT=$ROOT/tools
DATA_ROOT=$ROOT/data
PYTHON_ROOT=$TOOLS_ROOT/python
PERL_ROOT=$TOOLS_ROOT/perl
R_ROOT=$TOOLS_ROOT/r

function compile {
# $1 is the directory to compile in
    cd $1
    ./configure
    make
}

function download_gz {
    # $1 is the URL
    # $2 is the directory
    mkdir -p $2\
      && curl $1\
      | tar -xz --strip-components=1 -C $2
}

function download_zip {
    # $1 is the URL
    # $2 is the directory to extract into
    # $3 is the file name
    mkdir -p $2\
      && ZIP_FILE=$2/$3.zip\
      && wget $1 -O $ZIP_FILE\
      && unzip -d $2 $ZIP_FILE\
      && rm $ZIP_FILE
}

# Script Start

cd $TOOLS_ROOT

# Compile python and use it instead of the system python
download_gz https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz $PYTHON_ROOT

download_gz http://www.cpan.org/src/5.0/perl-$PERL_VERSION.tar.gz $PERL_ROOT\
  && mv $PERL_ROOT/configure.gnu $PERL_ROOT/configure.sh

download_gz http://cran.csiro.au/src/base/R-3/R-$R_VERSION.tar.gz $R_ROOT

download_zip https://dl.bintray.com/groovy/maven/apache-groovy-binary-$GROOVY_VERSION.zip $TOOLS_ROOT groovy\
  && mv $TOOLS_ROOT/groovy-$GROOVY_VERSION $TOOLS_ROOT/groovy

git clone https://github.com/ssadedin/groovy-ngs-utils --depth=1
