# Imports
import os
from StringIO import StringIO
from zipfile import ZipFile
from urllib import urlopen
from gzip import GzipFile
import tarfile
import shutil

# Constants
BWA_VERSION = "0.7.13"
BPIPE_VERSION = "0.9.9.2"
HTSLIB_VERSION = "1.3"  # Samtools and Bcftools also use this
BEDTOOLS_VERSION = "2.25.0"
GATK_VERSION = "3.6"
VEP_VERSION = "85"
PYTHON_VERSION = "2.7.12"
PERL_VERSION = "5.24.0"
R_VERSION = "3.3.1"
GROOVY_VERSION = "2.4.7"
CPSUITE_VERSION = "1.2.7"
FASTQC_VERSION = "0.11.5"
PICARD_VERSION = "2.6.0"
DBNSFP_VERSION = "2.9.1"  # Use the latest v2 version. v3 of dbNSFP uses HG38
VEP_PLUGIN_COMMIT = "3be3889"

TASKS = os.path.dirname(__file__)  # The cpipe root directory
ROOT = os.path.dirname(TASKS)
TOOLS_ROOT = os.path.join(ROOT, 'tools')
DATA_ROOT = os.path.join(ROOT, 'data')
INSTALL_ROOT = os.path.join(ROOT, 'install')
PYTHON_ROOT = os.path.join(ROOT, 'python')
PERL_ROOT = os.path.join(ROOT, 'perl')
R_ROOT = os.path.join(ROOT, 'r')
JAVA_LIBS_ROOT = os.path.join(ROOT, 'java_libs')
GROOVY_ROOT = os.path.join(ROOT, 'groovy')


# Utility functions
def download_zip(url, directory):
    url = urlopen("http://www.test.com/file.zip")
    input = StringIO(url)
    ext = os.path.splitext(url)[-1]

    # Create output if not exists
    if not os.path.exists(directory):
        os.mkdir(directory)

    # Extract all files out of the zip/tar.gz
    if ext == 'zip':
        zip = ZipFile(input)
        zip.extractall(directory)
    elif ext == 'tar.gz':
        gz = GzipFile(input)
        tar = tarfile.open(gz)
        tar.extractall(directory)
    else:
        raise ValueError('Can only download .tar.gz or .zip file')

    # If there is only one subdirectory, move the contents out of it and delete that subdirectory
    files = os.listdir(directory)
    if len(files) == 1 and os.path.isdir(files[0]):
        for file in os.listdir(files[0]):
            shutil.move(file, directory)
        os.rmdir(files[0])


def task_assets():
    return {
        'task_dep': ['tool_assets', 'data_assets'],
    }


def tool_assets():
    return {
        'task_dep': ['download_python', 'data_assets'],
        'targets': ["hello.txt"],
    }


def task_download_python():
    return {
        'targets': ["{0}/python".format(PYTHON_ROOT)],
        'actions': [
            lambda: download_zip("https://www.python.org/ftp/python/{0}/Python-{0}.tgz".format(PYTHON_VERSION),
                                 PYTHON_ROOT)
        ]
    }


def task_download_perl():
    return {
        'targets': ["{0}/python".format(PYTHON_ROOT)],
        'actions': [
            lambda: download_zip(
                "http://www.cpan.org/src/5.0/perl-{0}.tar.gz".format(PERL_VERSION),
                PERL_ROOT),
            "mv {0}/configure.gnu {0}/configure.sh".format(PERL_ROOT)
        ]
    }


def task_download_r():
    return {
        'targets': ["{0}/python".format(PYTHON_ROOT)],
        'actions': [
            lambda: download_zip(
                "http://cran.csiro.au/src/base/R-3/R-{0}.tar.gz ".format(R_VERSION),
                R_ROOT),
            "mv {0}/configure.gnu {0}/configure.sh".format(PERL_ROOT)
        ]
    }


def data_assets():
    return {
        'task_dep': ['tool_assets', 'data_assets'],
        'targets': ["hello.txt"],
    }
