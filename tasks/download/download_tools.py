'''
The download_* tasks are different from the install tasks because they
 - Do something that is platform independent (no compiling C code)
 - Remain all in one directory so we can zip it up
'''
from tasks.common import *
from tasks.nectar.nectar_util import *


def download_task(url, type=None):
    def action():
        temp_dir = tempfile.mkdtemp()
        download_zip(url, temp_dir, type=type)
        return {
            'dir': temp_dir
        }

    return {
        'actions': [action],
        'uptodate': [False]
    }


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
    return download_task(
        'http://apache.mirror.serversaustralia.com.au/maven/maven-3/{version}/binaries/apache-maven-{version}-bin.tar.gz'.format(
            version=MAVEN_VERSION))


def task_download_cpanm():
    if swift_install():
        return nectar_download('cpanm')
    else:
        def action():
            temp_dir = tempfile.mkdtemp()
            sh('''
                    curl -L https://cpanmin.us/ -o cpanm
                    chmod +x cpanm
            ''', cwd=temp_dir)

        return {
            'actions': [action],
            'uptodate': [False]
        }


def task_download_perl():
    if swift_install():
        return nectar_download('perl')
    else:
        return download_task("http://www.cpan.org/src/5.0/perl-{0}.tar.gz".format(PERL_VERSION))


def task_download_r():
    if swift_install():
        return nectar_download('r')
    else:
        return download_task("http://cran.csiro.au/src/base/R-3/R-{0}.tar.gz".format(R_VERSION))


def task_download_groovy():
    if swift_install():
        return nectar_download('groovy')
    else:
        return download_task("https://dl.bintray.com/groovy/maven/apache-groovy-binary-{0}.zip".format(GROOVY_VERSION), 'zip')


def task_download_bwa():
    if swift_install():
        return nectar_download('bwa')
    else:
        return download_task("https://codeload.github.com/lh3/bwa/tar.gz/v{0}".format(BWA_VERSION))


def task_download_htslib():
    if swift_install():
        return nectar_download('htslib')
    else:
        return download_task("https://github.com/samtools/htslib/releases/download/{0}/htslib-{0}.tar.bz2".format(
            HTSLIB_VERSION))


def task_download_samtools():
    if swift_install():
        return nectar_download('samtools')
    else:
        return download_task("https://github.com/samtools/samtools/releases/download/{0}/samtools-{0}.tar.bz2".format(
            HTSLIB_VERSION))


def task_download_bcftools():
    if swift_install():
        return nectar_download('bcftools')
    else:
        return download_task('https://github.com/samtools/bcftools/releases/download/{0}/bcftools-{0}.tar.bz2'.format(
            HTSLIB_VERSION))


def task_download_bedtools():
    if swift_install():
        return nectar_download('bedtools')
    else:
        return download_task("https://codeload.github.com/arq5x/bedtools2/tar.gz/v{0}".format(BEDTOOLS_VERSION))

def task_download_vep():
    if swift_install():
        return nectar_download('vep')
    else:
        return download_task("https://github.com/Ensembl/ensembl-vep/archive/release/{0}.tar.gz".format(VEP_VERSION))

def task_download_fastqc():
    if swift_install():
        return nectar_download('fastqc')
    else:
        return download_task(
            "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v{0}.zip".format(FASTQC_VERSION))

def task_download_bpipe():
    if swift_install():
        return nectar_download('bpipe')
    else:
        def action():
            temp_dir = tempfile.mkdtemp()
            sh('''
            git clone -c advice.detachedHead=false -b {bpipe_ver} --depth 1 https://github.com/ssadedin/bpipe {bpipe_dir}
            cd {bpipe_dir}
            ./gradlew dist
            '''.format(bpipe_ver=BPIPE_VERSION, bpipe_dir=temp_dir))
            return {'dir': temp_dir}

        return {
            'actions': [action],
            'uptodate': [False]
        }


def task_download_gatk():
    if swift_install():
        return nectar_download('gatk')
    else:
        def action():
            temp_dir = tempfile.mkdtemp()
            download_zip("https://codeload.github.com/broadgsa/gatk-protected/tar.gz/{}".format(GATK_VERSION),
                         temp_dir,
                         type='tgz')
            sh('''
                mvn verify -P\!queue
                GATK_JAR=`readlink -f target/GenomeAnalysisTK.jar`
                unlink target/GenomeAnalysisTK.jar
                mv $GATK_JAR ./GenomeAnalysisTK.jar
                bash -O extglob -O dotglob -c 'rm -rf !(GenomeAnalysisTK.jar)'
           ''', cwd=temp_dir)
            return {'dir': temp_dir}

        return {
            'actions': [action],
            'task_dep': ['install_maven'],
            'uptodate': [False]
        }


def task_download_picard():
    if swift_install():
       return nectar_download('picard')
    else:
        def action():
            temp_dir = tempfile.mkdtemp()
            urlretrieve(
                'https://github.com/broadinstitute/picard/releases/download/{0}/picard.jar'.format(PICARD_VERSION),
                os.path.join(temp_dir, 'picard.jar')
            )
            return {'dir': temp_dir}

        return {
            'actions': [action],
            'uptodate': [False]
        }


def task_download_perl_libs():
    """
    Downloads all cpan libs into the cpan directory. Each dependency is a subtask
    :return:
    """
    if swift_install():
        return nectar_download('perl_libs')
    else:
        def action():
            temp_dir = tempfile.mkdtemp()
            sh('''
            # Module::Build has to be installed to even work out the dependencies of perl modules, so we do that first
            # (while also saving the archive so Module::Build will be bundled on NECTAR)
            cpanm -l {perl_lib} --save-dists {cpan} Module::Build
            # Now, download archives of everything we need without installing them
            cpanm --save-dists {cpan} -L /dev/null --scandeps --installdeps .
            '''.format(cpan=temp_dir, perl_lib=PERL_LIB_ROOT))
            return {'dir': temp_dir}

        return {
            'task_dep': ['copy_config', 'download_perl', 'install_perl', 'download_cpanm'],
            'uptodate': [False],
            'actions': [action]
        }


def task_download_vep_libs():
    if swift_install():
        return nectar_download('vep_libs')
    else:
        def action():
            temp_dir = tempfile.mkdtemp()
            sh('perl {vep_dir}/INSTALL.pl --NO_HTSLIB --NO_TEST --AUTO a --DESTDIR {vep_libs}'.format(vep_dir=VEP_ROOT,
                                                                                            vep_libs=temp_dir))
            return {'dir': temp_dir}
        return {
            'task_dep': ['copy_config', 'install_perl_libs', 'install_perl', 'install_vep'],
            'actions': [action],
            'uptodate': [False]
        }


def task_download_vep_plugins():
    if swift_install():
        return nectar_download('vep_plugins')
    else:
        def action():
            temp_dir = tempfile.mkdtemp()
            sh('''
                git init
                git remote add origin https://github.com/Ensembl/VEP_plugins
                git fetch
                git checkout -t origin/master
                git reset --hard {vep_plugin_commit}
                rm -rf .git
            '''.format(vep_plugin_commit=VEP_PLUGIN_COMMIT), cwd=temp_dir)
            return {'dir': temp_dir}
        return {
            'actions': [action],
            'task_dep': ['copy_config'],
            'uptodate': [False]
        }


def task_download_java_libs():
    if swift_install():
        return nectar_download('java_libs')
    else:
        def action():
            temp_dir = tempfile.mkdtemp()
            sh('gradle copyDeps -Pdir={}'.format(temp_dir))
            return {'dir': temp_dir}
            
        return {
            'actions': [action],
            'uptodate': [False]
        }


def task_download_c_libs():
    return {
        'actions': None,
        'task_dep': ['download_bzip2', 'download_xz', 'download_pcre', 'download_libcurl', 'download_zlib']
    }


def task_download_bzip2():
    if swift_install():
        return nectar_download('bzip2')
    else:
        return download_task(
            'http://www.bzip.org/{0}/bzip2-{0}.tar.gz'.format(BZIP_VERSION))


def task_download_xz():
    if swift_install():
        return nectar_download('xz')
    else:
        return download_task(
            'http://tukaani.org/xz/xz-{}.tar.gz'.format(XZ_VERSION))


def task_download_pcre():
    if swift_install():
        return nectar_download('pcre')
    else:
        return download_task(
            'ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-{}.tar.gz'.format(PCRE_VERSION))


def task_download_libcurl():
    if swift_install():
        return nectar_download('libcurl')
    else:
        return download_task(
            'https://curl.haxx.se/download/curl-{}.tar.gz'.format(LIBCURL_VERSION))


def task_download_zlib():
    if swift_install():
        return nectar_download('zlib')
    else:
        return download_task(
            'https://codeload.github.com/madler/zlib/tar.gz/v{}'.format(ZLIB_VERSION))
