# Imports
import os
from StringIO import StringIO
from zipfile import ZipFile
from urllib import urlopen, urlretrieve
from ftplib import FTP
import tarfile
import shutil
import tempfile
import glob

from tasks.common import cmd, TOOLS_ROOT, DATA_ROOT, INSTALL_ROOT
from tasks.manual_assets.dependencies import *

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

PYTHON_ROOT = os.path.join(TOOLS_ROOT, 'python')
PERL_ROOT = os.path.join(TOOLS_ROOT, 'perl')
R_ROOT = os.path.join(TOOLS_ROOT, 'r')
JAVA_LIBS_ROOT = os.path.join(TOOLS_ROOT, 'java_libs')
GROOVY_ROOT = os.path.join(TOOLS_ROOT, 'groovy')

VEP_ROOT = os.path.join(TOOLS_ROOT, 'vep')
PERL_LIB_ROOT = os.path.join(TOOLS_ROOT, 'perl_lib')


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


def download_ftp_list(ftp, files, target_dir):
    ftp.login()
    for file in files:
        ftp.retrbinary(
            'RETR {}'.format(file),
            open(os.path.join(target_dir, file), 'wb').write
        )


def task_manual_assets():
    return {
        'actions': None,
        'task_dep': ['manual_install_dependencies', 'tool_assets', 'data_assets']
    }


def task_tool_assets():
    return {
        'actions': None,
        'task_dep': [
            'download_python',
            'download_perl',
            'download_r',
            'download_groovy',
            'download_bwa',
            'download_htslib',
            'download_samtools',
            'download_bcftools',
            'download_bedtools',
            'download_vep',
            'download_fastqc',
            'download_bpipe',
            'download_gatk',
            'download_picard',
            'download_perl_libs',
            'download_vep_libs',
            'download_vep_plugins',
            'download_java_libs'
        ],
    }


def task_download_python():
    return {
        'targets': [PYTHON_ROOT],
        'actions': [
            lambda: download_zip("https://www.python.org/ftp/python/{0}/Python-{0}.tgz".format(PYTHON_VERSION),
                                 PYTHON_ROOT)
        ],
        'uptodate': [True]
    }


def task_download_perl():
    return {
        'targets': [PERL_ROOT],
        'actions': [
            lambda: download_zip("http://www.cpan.org/src/5.0/perl-{0}.tar.gz".format(PERL_VERSION), PERL_ROOT),
            "mv {0}/configure.gnu {0}/configure.sh".format(PERL_ROOT),
            lambda: os.environ['']
        ],
        'uptodate': [True]
    }


def task_download_r():
    return {
        'targets': [R_ROOT],
        'actions': [
            lambda: download_zip("http://cran.csiro.au/src/base/R-3/R-{0}.tar.gz".format(R_VERSION), R_ROOT),
        ],
        'uptodate': [True]
    }


def task_download_groovy():
    return {
        'targets': [GROOVY_ROOT],
        'actions': [
            lambda: download_zip(
                "https://dl.bintray.com/groovy/maven/apache-groovy-binary-{0}.zip".format(GROOVY_VERSION), GROOVY_ROOT)
        ],
        'uptodate': [True]
    }


def task_download_bwa():
    BWA_ROOT = os.path.join(TOOLS_ROOT, 'bwa')
    return {
        'targets': [BWA_ROOT],
        'actions': [
            lambda: download_zip(
                "https://codeload.github.com/lh3/bwa/tar.gz/v{0}".format(BWA_VERSION), BWA_ROOT, type='tgz')
        ],
        'uptodate': [True]
    }


def task_download_htslib():
    HTSLIB_ROOT = os.path.join(TOOLS_ROOT, 'htslib')
    return {
        'targets': [HTSLIB_ROOT],
        'actions': [
            lambda: download_zip(
                "https://codeload.github.com/samtools/htslib/tar.gz/{0}".format(HTSLIB_VERSION), HTSLIB_ROOT,
                type='tgz')
        ],
        'uptodate': [True]
    }


def task_download_samtools():
    SAMTOOLS_ROOT = os.path.join(TOOLS_ROOT, 'samtools')
    return {
        'targets': [SAMTOOLS_ROOT],
        'actions': [
            lambda: download_zip(
                "https://codeload.github.com/samtools/samtools/tar.gz/{0}".format(HTSLIB_VERSION), SAMTOOLS_ROOT,
                type='tgz')
        ],
        'uptodate': [True]
    }


def task_download_bcftools():
    BCFTOOLS_ROOT = os.path.join(TOOLS_ROOT, 'bcftools')
    return {
        'targets': [BCFTOOLS_ROOT],
        'actions': [
            lambda: download_zip(
                "https://codeload.github.com/samtools/bcftools/tar.gz/{0}".format(HTSLIB_VERSION), BCFTOOLS_ROOT,
                type='tgz')
        ],
        'uptodate': [True]
    }


def task_download_bedtools():
    BEDTOOLS_ROOT = os.path.join(TOOLS_ROOT, 'bedtools')
    return {
        'targets': [BEDTOOLS_ROOT],
        'actions': [
            lambda: download_zip(
                "https://codeload.github.com/arq5x/bedtools2/tar.gz/v{0}".format(BEDTOOLS_VERSION), BEDTOOLS_ROOT,
                type='tgz')
        ],
        'uptodate': [True]
    }


def task_download_vep():
    def extract_vep():
        """Download ensembl tools and move VEP out of the tools directory"""

        VEP_TEMP = tempfile.mkdtemp()
        VEP_SUBDIR = os.path.join(VEP_TEMP, 'scripts', 'variant_effect_predictor')

        download_zip("https://github.com/Ensembl/ensembl-tools/archive/release/{0}.zip".format(VEP_VERSION), VEP_TEMP)
        os.rename(VEP_SUBDIR, VEP_ROOT)
        shutil.rmtree(VEP_TEMP)

    return {
        'targets': [VEP_ROOT],
        'actions': [extract_vep],
        'uptodate': [True]
    }


def task_download_fastqc():
    FASTQC_ROOT = os.path.join(TOOLS_ROOT, 'fastqc')
    return {
        'targets': [FASTQC_ROOT],
        'actions': [
            lambda: download_zip(
                "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v{0}.zip".format(FASTQC_VERSION),
                FASTQC_ROOT)
        ],
        'uptodate': [True]
    }


def task_download_bpipe():
    BPIPE_ROOT = os.path.join(TOOLS_ROOT, 'bpipe')
    return {
        'targets': [BPIPE_ROOT],
        'actions': [cmd('''
            git clone -b {bpipe_ver} --depth 1 https://github.com/ssadedin/bpipe {bpipe_dir}\
            && cd {bpipe_dir}\
            && ./gradlew dist
            '''.format(bpipe_dir=BPIPE_ROOT, bpipe_ver=BPIPE_VERSION), cwd=TOOLS_ROOT
                        )],
        'uptodate': [True]
    }


def task_download_gatk():
    GATK_ROOT = os.path.join(TOOLS_ROOT, 'gatk')
    return {
        'targets': [GATK_ROOT],
        'actions': [
            lambda: download_zip("https://codeload.github.com/broadgsa/gatk-protected/tar.gz/{}".format(GATK_VERSION),
                                 GATK_ROOT,
                                 type='tgz')
        ],
        'uptodate': [True]
    }


def task_download_picard():
    PICARD_ROOT = os.path.join(TOOLS_ROOT, 'picard')

    def action():
        os.makedirs(PICARD_ROOT)
        urlretrieve(
            'https://github.com/broadinstitute/picard/releases/download/{0}/picard.jar'.format(PICARD_VERSION),
            os.path.join(PICARD_ROOT, 'picard.jar')
        )

    return {
        'targets': [os.path.join(PICARD_ROOT, 'picard.jar')],
        'actions': [action],
        'uptodate': [True]
    }


def task_download_perl_libs():
    return {
        'targets': [PERL_LIB_ROOT],
        'actions': [
            lambda: os.makedirs(PERL_LIB_ROOT),
            cmd('cpanm --installdeps --local-lib-contained {}/perl_lib .'.format(TOOLS_ROOT), cwd=INSTALL_ROOT)
        ],
        'uptodate': [True]
    }


def task_download_vep_libs():
    return {
        'targets': [os.path.join(PERL_LIB_ROOT, 'lib/perl5/Bio/EnsEMBL')],
        'actions': [
            cmd('yes | perl {vep_dir}/INSTALL.pl --NO_HTSLIB --AUTO a --DESTDIR {perl5}'.format(
                vep_dir=VEP_ROOT,
                perl5=os.path.join(PERL_LIB_ROOT, 'lib/perl5')
            ))
        ],
        'uptodate': [True]
    }


def task_download_vep_plugins():
    VEP_PLUGIN_ROOT = os.path.join(TOOLS_ROOT, 'vep_plugins')
    return {
        'targets': [os.path.join(VEP_PLUGIN_ROOT, 'Condel.pm')],
        'actions': [
            cmd('''
            git init\
            && git remote add origin https://github.com/Ensembl/VEP_plugins\
            && git fetch\
            && git checkout -t origin/master\
            && git reset --hard $VEP_PLUGIN_COMMIT\
            && rm -rf .git
            ''', cwd=VEP_PLUGIN_ROOT)
        ],
        'uptodate': [True]
    }


def task_download_java_libs():
    return {
        'targets': [JAVA_LIBS_ROOT],
        'actions': None,
        'task_dep': [
            'download_junit_xml_formatter',
            'download_groovy_ngs_utils',
            'download_takari_cpisuite'
        ]
    }


def task_download_junit_xml_formatter():
    return {
        'actions': [
            cmd('''
                git clone https://github.com/barrypitman/JUnitXmlFormatter\
                && pushd JUnitXmlFormatter\
                    && mvn install\
                    && mv target/JUnitXmlFormatter* {java_libs_dir}\
                && popd\
                && rm -rf JUnitXmlFormatter
            '''.format(java_libs_dir=JAVA_LIBS_ROOT), cwd=JAVA_LIBS_ROOT)
        ],
        'uptodate': [
            lambda: len(glob.glob(os.path.join(JAVA_LIBS_ROOT, 'JUnitXmlFormatter*.jar'))) > 0
        ]
    }


def task_download_groovy_ngs_utils():
    return {
        'targets': [os.path.join(JAVA_LIBS_ROOT, 'groovy-ngs-utils.jar')],
        'actions': [
            cmd('''
              git clone https://github.com/ssadedin/groovy-ngs-utils -b upgrade-biojava --depth=1 --quiet\
                && pushd groovy-ngs-utils\
                && ./gradlew jar \
                && popd\
                && mv {java_libs_dir}/groovy-ngs-utils/build/libs/groovy-ngs-utils.jar {java_libs_dir}\
                && rm -rf groovy-ngs-utils
            '''.format(java_libs_dir=JAVA_LIBS_ROOT), cwd=JAVA_LIBS_ROOT)
        ],
        'uptodate': [True],
    }


def task_download_takari_cpisuite():
    return {
        'actions': [
            cmd('''
             mvn dependency:copy \
                    -Dartifact=io.takari.junit:takari-cpsuite:{cpsuite_version}\
                    -DoutputDirectory={java_libs_dir}\
                    -DstripVersion=true
            '''.format(cpsuite_version=CPSUITE_VERSION, java_libs_dir=JAVA_LIBS_ROOT), cwd=JAVA_LIBS_ROOT)
        ],
        'uptodate': [
            lambda: len(glob.glob(os.path.join(JAVA_LIBS_ROOT, 'takari-cpsuite*'))) > 0
        ],
    }


def task_data_assets():
    return {
        'actions': None,
        'task_dep': [
            'download_junit_xml_formatter',
            'download_groovy_ngs_utils',
            'download_takari_cpisuite'
        ]
    }


def task_download_dbnsfp():
    DBNSFP_ROOT = os.path.join(DATA_ROOT, 'dbnsfp')

    return {
        'targets': [os.path.join(DATA_ROOT, 'dbnsfp', 'dbNSFP.gz')],
        'actions': [
            lambda: download_zip("ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv{}.zip".format(DBNSFP_VERSION),
                                 DBNSFP_ROOT),
            cmd('''
            source {install_dir}/environment.sh\
            && mkdir -p dbnsfp\
            && pushd {data_dir}/dbnsfp\
                && head -n1 dbNSFP*chr1 > h\
                && cat dbNSFP*chr* | grep -v ^#chr | sort -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP.gz\
                && tabix -s 1 -b 2 -e 2 dbNSFP.gz\
                && bash -O extglob -c 'rm -rf !(dbNSFP.gz*)'
            '''.format(install_dir=INSTALL_ROOT, data_dir=DATA_ROOT), cwd=DATA_ROOT, executable='bash')
        ],
        'task_dep': [
            'download_htslib'
        ],
        'uptodate': [True],
    }


def task_install_vep_cache():
    VEP_CACHE = os.path.join(DATA_ROOT, 'vep_cache')
    return {
        'targets': [VEP_CACHE],
        'actions': [
            lambda: os.makedirs(VEP_CACHE),
            '''perl {tools_dir}/vep/INSTALL.pl\
            --NO_HTSLIB\
            --CACHEDIR $VEP_CACHE\
            --AUTO cf\
            --SPECIES homo_sapiens_refseq\
            --ASSEMBLY GRCh37'''.format(tools_dir=TOOLS_ROOT)
        ],
        'task_dep': [
            'download_htslib'
        ],
        'uptodate': [True],
    }


def task_download_ucsc():
    UCSC_ROOT = os.path.join(DATA_ROOT, 'ucsc')

    return {
        'targets': [UCSC_ROOT],
        'actions': [
            lambda: download_ftp_list(
                FTP("ftp://ftp.broadinstitute.org/bundle/2.8/hg19/",
                    user="gsapubftp-anonymous:cpipe.user@cpipeline.org"),
                ["ucsc.hg19.dict.gz", "ucsc.hg19.fasta.gz", "ucsc.hg19.fasta.fai.gz"],
                UCSC_ROOT
            )
        ],
        'uptodate': [True],
    }


def task_download_mills_and_1000g():
    MILLS_ROOT = os.path.join(DATA_ROOT, 'mills_and_1000g')

    return {
        'targets': [MILLS_ROOT],
        'actions': [
            lambda: download_ftp_list(
                FTP("ftp://ftp.broadinstitute.org/bundle/2.8/hg19/",
                    user="gsapubftp-anonymous:cpipe.user@cpipeline.org"),
                ["Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
                 "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz"],
                MILLS_ROOT
            )
        ],
        'uptodate': [True],
    }


def task_download_dbsnp():
    DBSNP_ROOT = os.path.join(DATA_ROOT, 'dbsnp')

    return {
        'targets': [DBSNP_ROOT],
        'actions': [
            lambda: download_ftp_list(
                FTP("ftp://ftp.broadinstitute.org/bundle/2.8/hg19/",
                    user="gsapubftp-anonymous:cpipe.user@cpipeline.org"),
                ["dbsnp_138.hg19.vcf.gz", "dbsnp_138.hg19.vcf.idx.gz"],
                DBSNP_ROOT
            )
        ],
        'uptodate': [True],
    }


TRIO_REFINEMENT_FILE = '{data_dir}/1000G_phase3/1000G_phase3_v4_20130502.sites.hg19.vcf.gz'.format(data_dir=DATA_ROOT)


def task_convert_trio_refinement():
    return {
        'actions': [
            cmd(
                '''
                 source {install_dir}/environment.sh
                 java -jar $TOOLS_ROOT/picard/picard.jar LiftoverVcf \
                        I={data_dir}/1000G_phase3/1000G_phase3_v4_20130502.sites.vcf \
                        O={data_dir}/1000G_phase3/1000G_phase3_v4_20130502.sites.hg19.vcf \
                        CHAIN={data_dir}/1000G_phase3/b37tohg19.chain \
                        REJECT={data_dir}/1000G_phase3/liftover.rejected_variants.vcf \
                        R={data_dir}/ucsc/ucsc.hg19.fasta\
                    && bgzip {data_dir}/1000G_phase3/1000G_phase3_v4_20130502.sites.hg19.vcf \
                    && tabix -p vcf {data_dir}/1000G_phase3/1000G_phase3_v4_20130502.sites.hg19.vcf.gz

                bash -O extglob -c "rm -rf {data_dir}/1000G_phase3/!(1000G_phase3_v4_20130502.sites.hg19.vcf.gz*)"
                '''.format(data_dir=DATA_ROOT, install_dir=INSTALL_ROOT),
                executable='bash'
            ),
            ''.format(DATA_ROOT)
        ],
        'targets': [TRIO_REFINEMENT_FILE],
        'task_dep': [
            'download_trio_refinement',
            'download_refinement_liftover',
        ],
        'uptodate': [True]
    }


def task_download_trio_refinement():
    RAW_VCF = '{}/1000G_phase3/1000G_phase3_v4_20130502.sites.vcf'.format(DATA_ROOT)

    return {
        'actions': [
            '''
            mkdir -p {data_dir}/1000G_phase3\
            && curl \
                --user gsapubftp-anonymous:cpipe.user@cpipeline.org \
                 ftp://ftp.broadinstitute.org/bundle/2.8/b37/1000G_phase3_v4_20130502.sites.vcf.gz \
                 | gunzip > {data_dir}/1000G_phase3/1000G_phase3_v4_20130502.sites.vcf
            '''.format(data_dir=DATA_ROOT)
        ],
        'targets': [RAW_VCF],
        # The task is up to date if the final refinement file exists or if just this step's product exists
        'uptodate': [lambda: os.path.exists(TRIO_REFINEMENT_FILE) or os.path.exists(RAW_VCF)]
    }


def task_download_refinement_liftover():
    LIFTOVER_FILE = '{}/b37tohg19.chain'.format(DATA_ROOT)

    return {
        'actions': [
            '''
           wget \
                    --user gsapubftp-anonymous:cpipe.user@cpipeline.org \
                     -P {data_dir}/1000G_phase3 \
                     ftp://gsapubftp-anonymous@ftp.broadinstitute.org/Liftover_Chain_Files/b37tohg19.chain
            '''.format(data_dir=DATA_ROOT)
        ],
        'targets': [LIFTOVER_FILE],
        # The task is up to date if the final refinement file exists or if just this step's product exists
        'uptodate': [lambda: os.path.exists(TRIO_REFINEMENT_FILE) or os.path.exists(LIFTOVER_FILE)]
    }


def task_download_chromosome_sizes():
    CHROMOSOME_DIR = os.path.join(DATA_ROOT, 'chromosomes')
    CHROMOSOME_FILE = os.path.join(CHROMOSOME_DIR, 'hg19.genome')

    return {
        'targets': [CHROMOSOME_FILE],
        'actions': [
            lambda: os.makedirs(CHROMOSOME_DIR),
            '''mysql
            --user=genome\
            --host=genome-mysql.cse.ucsc.edu\
            -A\
            -e "select chrom, size from hg19.chromInfo"\
            > {data}/chromosomes/hg19.genome'''.format(data=DATA_ROOT)
        ]
    }


def task_index_reference_files():
    return {
        'task_dep': [
            'bwa_index_ucsc_reference',
            'samtools_index_ucsc_reference'
        ],
        'actions': None
    }


def task_bwa_index_ucsc_reference():
    UCSC_BWA_INDEX = '{data}/ucsc/ucsc.hg19.fasta.bwt'.format(data=DATA_ROOT)

    return {
        'targets': [UCSC_BWA_INDEX],
        'actions': [
            '{tools}/bwa/bwa index -a bwtsw {data}/ucsc/ucsc.hg19.fasta'.format(tools=TOOLS_ROOT, data=DATA_ROOT)
        ],
        'task_dep': [
            'download_bwa',
            'download_ucsc'
        ],
        'uptodate': [True]
    }


def task_samtools_index_ucsc_reference():
    UCSC_SAMTOOLS_INDEX = '{data}/ucsc/ucsc.hg19.fasta.fai'.format(data=DATA_ROOT)

    return {
        'targets': [UCSC_SAMTOOLS_INDEX],
        'actions': [
            '{tools}/samtools/samtools faidx {data}/ucsc/ucsc.hg19.fasta'.format(tools=TOOLS_ROOT, data=DATA_ROOT)
        ],
        'task_dep': [
            'download_samtools',
            'download_ucsc'
        ],
        'uptodate': [True]
    }
