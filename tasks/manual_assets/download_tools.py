from StringIO import StringIO
from zipfile import ZipFile
from urllib import urlopen, urlretrieve
import tarfile
import shutil
import tempfile
import glob

from tasks.common import *

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

