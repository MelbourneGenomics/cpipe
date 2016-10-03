import os
import tarfile
import shutil
from urllib import urlopen, urlretrieve
from StringIO import StringIO
from doit.action import CmdAction
from zipfile import ZipFile

# General paths
HERE = os.path.dirname(__file__)  # The cpipe root directory
ROOT = os.path.dirname(HERE)
TOOLS_ROOT = os.path.join(ROOT, 'tools')
DATA_ROOT = os.path.join(ROOT, 'data')
INSTALL_ROOT = os.path.join(ROOT, 'install')

# Versions
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

# Tool paths
PYTHON_ROOT = os.path.join(TOOLS_ROOT, 'python')
PERL_ROOT = os.path.join(TOOLS_ROOT, 'perl')
R_ROOT = os.path.join(TOOLS_ROOT, 'r')
JAVA_LIBS_ROOT = os.path.join(TOOLS_ROOT, 'java_libs')
GROOVY_ROOT = os.path.join(TOOLS_ROOT, 'groovy')
VEP_ROOT = os.path.join(TOOLS_ROOT, 'vep')
PERL_LIB_ROOT = os.path.join(TOOLS_ROOT, 'perl_lib')
BWA_ROOT = os.path.join(TOOLS_ROOT, 'bwa')
HTSLIB_ROOT = os.path.join(TOOLS_ROOT, 'htslib')
SAMTOOLS_ROOT = os.path.join(TOOLS_ROOT, 'samtools')
BCFTOOLS_ROOT = os.path.join(TOOLS_ROOT, 'bcftools')
BEDTOOLS_ROOT = os.path.join(TOOLS_ROOT, 'bedtools')
GATK_ROOT = os.path.join(TOOLS_ROOT, 'gatk')

ENVIRONMENT_FILE = os.path.join(INSTALL_ROOT, 'environment.sh')

# Utility functions
def download_zip(url_str, directory, type=None):
    url = urlopen(url_str)
    input = StringIO(url.read())

    # Create output if not exists
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Either take the type from the arguments, or else deduce it from the URL
    if not type:
        [name, ext2] = os.path.splitext(url_str)
        [name, ext1] = os.path.splitext(name)

        if ext2 == '.zip':
            type = 'zip'
        elif (ext1 == '.tar' and ext2 == '.gz') or ext2 == '.tgz':
            type = 'tgz'

    if type == 'zip':
        zip = ZipFile(input)
        zip.extractall(directory)
    elif type == 'tgz':
        tar = tarfile.open(fileobj=input, mode='r:gz')
        tar.extractall(directory)
    else:
        raise ValueError('Can only download .tar.gz or .zip file')

    # If there is only one subdirectory, move the contents out of it and delete that subdirectory
    files = [os.path.join(directory, f) for f in os.listdir(directory)]
    if len(files) == 1 and os.path.isdir(files[0]):
        subdir = files[0]
        for subfile in [os.path.join(subdir, f) for f in os.listdir(subdir)]:
            shutil.move(subfile, directory)
        os.rmdir(files[0])


def cmd(command, **kwargs):
    """
    Creates a doit CmdAction with certain default parameters used by cpipe, namely using bash as the default shell,
    using the cpipe root as the working directory, and sourcing the environment file before use. Using the environment
    file instead of simply setting the subprocess's env is done to prevent the environment file becoming out of sync
    with the actual environment variables we need set
    :param command: The bash command to run
    :param kwargs: And additional arguments to pass to CmdAction or popen
    :return: A doit CmdAction
    """
    return CmdAction(
        'source {}\n'.format(ENVIRONMENT_FILE) + command,
        executable='bash',
        shell=True,
        cwd=ROOT,
        **kwargs
    )