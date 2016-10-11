from doit.tools import create_folder
import tempfile
import glob
import re
from urllib import urlopen, urlretrieve

from tasks.common import *


def task_tool_assets():
    """
    Downloads all the tools needed to run cpipe and install other tools
    :return:
    """
    return {
        'actions': None,
        'task_dep': [
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
            'download_java_libs',
        ],
    }


def task_download_maven():
    return {
        'targets': [MAVEN_ROOT],
        'actions': [
            lambda: create_folder(MAVEN_ROOT),
            lambda: download_zip(
                'http://apache.mirror.serversaustralia.com.au/maven/maven-3/{version}/binaries/apache-maven-{version}-bin.tar.gz'.format(
                    version=MAVEN_VERSION),
                MAVEN_ROOT
            )
        ],
        'uptodate': [True]
    }


def task_download_cpanm():
    return {
        'targets': [CPANM_ROOT],
        'actions': [
            lambda: create_folder(CPANM_ROOT),
            cmd('''
                curl -L https://cpanmin.us/ -o cpanm
                chmod +x cpanm
            ''', cwd=CPANM_ROOT),
        ],
        'task_dep': ['copy_config'],
        'uptodate': [True]
    }


def task_download_perl():
    return {
        'targets': [PERL_ROOT],
        'actions': [
            lambda: download_zip("http://www.cpan.org/src/5.0/perl-{0}.tar.gz".format(PERL_VERSION), PERL_ROOT),
            "mv {0}/configure.gnu {0}/configure.sh".format(PERL_ROOT),
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
    FASTQC_EXE = os.path.join(FASTQC_ROOT, 'fastqc')

    return {
        'targets': [FASTQC_ROOT],
        'actions': [
            lambda: download_zip(
                "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v{0}.zip".format(FASTQC_VERSION),
                FASTQC_ROOT),
            'chmod +x {fastqc}'.format(fastqc=FASTQC_EXE)
        ],
        'uptodate': [True]
    }


def task_download_bpipe():
    BPIPE_ROOT = os.path.join(TOOLS_ROOT, 'bpipe')
    return {
        'targets': [BPIPE_ROOT],
        'actions': [
            cmd('''
            git clone -c advice.detachedHead=false -b {bpipe_ver} --depth 1 https://github.com/ssadedin/bpipe {bpipe_dir}\
            && cd {bpipe_dir}\
            && ./gradlew dist
            '''.format(bpipe_dir=BPIPE_ROOT, bpipe_ver=BPIPE_VERSION), cwd=TOOLS_ROOT)
        ],
        'task_dep': ['copy_config'],
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
        create_folder(PICARD_ROOT)
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
    """
    Downloads all cpan libs into the cpan directory. Each dependency is a subtask
    :return:
    """

    return {
        'targets': [CPAN_ROOT],
        'task_dep': ['copy_config', 'download_perl', 'compile_perl', 'download_cpanm'],
        'uptodate': [True],
        'actions': [
            lambda: create_folder(CPAN_ROOT),
            lambda: create_folder(PERL_LIB_ROOT),

            # Module::Build has to be installed to even work out the dependencies of perl modules, so we do that first
            # (while also saving the archive so Module::Build will be bundled on NECTAR)
            cmd('cpanm -l {perl_lib} --save-dists {cpan} Module::Build'.format(perl_lib=PERL_LIB_ROOT, cpan=CPAN_ROOT),
                env=get_cpanm_env()),

            # Now, download archives of everything we need without installing them
            cmd('cpanm --save-dists {cpan} -L /dev/null --scandeps --installdeps .'.format(cpan=CPAN_ROOT),
                cwd=ROOT, env=get_cpanm_env())
        ]
    }


def task_download_vep_libs():
    return {
        'targets': [VEP_LIBS_ROOT],
        'task_dep': ['copy_config', 'install_perl_libs'],
        'actions': [
            cmd('yes | perl {vep_dir}/INSTALL.pl --NO_HTSLIB --AUTO a --DESTDIR {vep_libs}'.format(
                vep_dir=VEP_ROOT,
                vep_libs=VEP_LIBS_ROOT
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
        'task_dep': ['copy_config'],
        'uptodate': [True]
    }


def task_download_java_libs():
    return {
        'actions': None,
        'task_dep': [
            'download_junit_xml_formatter',
            'download_groovy_ngs_utils',
            'download_takari_cpisuite'
        ]
    }


def task_make_java_libs_dir():
    return {
        'actions': [
            create_folder(JAVA_LIBS_ROOT),
        ],
        'targets': [JAVA_LIBS_ROOT],
        'uptodate': [True]
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
        'task_dep': ['copy_config', 'make_java_libs_dir', 'download_maven'],
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
        'task_dep': ['copy_config', 'make_java_libs_dir', 'download_maven'],
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
        'task_dep': ['copy_config', 'make_java_libs_dir', 'download_maven'],
        'uptodate': [
            lambda: len(glob.glob(os.path.join(JAVA_LIBS_ROOT, 'takari-cpsuite*'))) > 0
        ],
    }


def task_download_bzip2():
    def action():
        create_folder(BZIP_ROOT)
        download_zip(
            'http://www.bzip.org/{0}/bzip2-{0}.tar.gz'.format(BZIP_VERSION),
            BZIP_ROOT
        )

    return {
        'targets': [BZIP_ROOT],
        'actions': [action],
        'uptodate': [True]
    }

def task_download_xz():
    def action():
        create_folder(XZ_ROOT)
        download_zip(
            'http://tukaani.org/xz/xz-{}.tar.gz'.format(XZ_VERSION),
            XZ_ROOT
        )

    return {
        'targets': [XZ_ROOT],
        'actions': [action],
        'uptodate': [True]
    }

