from tasks.common import *
from doit.tools import create_folder
from tasks.compile_c_libs import *
import os
import subprocess


def task_install_all():
    'Compiles all assets including java assets; this is only needed for a manual install'
    return {
        'actions': None,
        'task_dep': [
            'install_perl',
            'install_r',
            'install_bwa',
            'install_htslib',
            'install_samtools',
            'install_bcftools',
            'install_bedtools',
            'install_gatk',
            'install_perl_libs'
        ]
    }


def task_install_nectar():
    'Compiles assets downloaded from nectar; excludes java compilation since this should already have been done'
    return {
        'actions': None,
        'task_dep': [
            'install_perl',
            'install_r',
            'install_bwa', 'install_htslib',
            'install_samtools',
            'install_bcftools',
            'install_bedtools',
            'install_perl_libs'
        ]
    }


def task_install_perl():
    def action(perl_root):
        sh('''
              ./Configure -de -Dprefix={}
              make
              make install
        '''.format(INSTALL_ROOT),
           cwd=perl_root
           )

    return {
        'actions': [action],
        'getargs': {'perl_dir': ('download_perl', 'dir')},
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_perl'],
        'targets': [os.path.join(INSTALL_ROOT, 'bin', 'perl')],
        'uptodate': [True]
    }


def task_install_r():
    def action(r_dir):
        sh('''
            ./configure
             make
             make prefix={} install
        '''.format(r_dir), cwd=r_dir)

    return {
        'actions': [action],
        'task_dep': ['download_perl', 'download_r', 'install_bzip2', 'install_xz', 'install_pcre', 'install_libcurl',
                     'install_zlib'],
        'getargs': {'r_dir': ('download_r', 'dir')},
        'targets': [os.path.join(INSTALL_ROOT, 'bin', 'R')],
        'uptodate': [True]
    }


def task_install_bwa():
    def action(bwa_dir):
        sh('make', cwd=bwa_dir)

    return {
        'actions': [action],
        'getargs': {'bwa_dir': ('download_bwa', 'dir')},
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_bwa'],
        'targets': [os.path.join(INSTALL_BIN, 'bwa')],
        'uptodate': [True]
    }


def task_install_htslib():
    def action(htslib_dir):
        sh('''
        ./configure --prefix={0}
        make
        make prefix={0} install
        '''.format(INSTALL_ROOT), cwd=htslib_dir)

    return {
        'task_dep': ['download_htslib', 'install_zlib'],
        'actions': [action],
        'getargs': {'htslib_dir': ('download_htslib', 'dir')},
        'targets': [os.path.join(INSTALL_ROOT, 'bin', 'htsfile')],
        'uptodate': [True]
    }


def task_install_samtools():
    def action(samtools_dir):
        sh('''
                ./configure --prefix={0} --with-htslib={1}
                make
                make prefix={0} install
        '''.format(INSTALL_ROOT, INSTALL_LIB), cwd=samtools_dir)

    return {
        'task_dep': ['install_zlib', 'install_htslib', 'download_samtools'],
        'actions': [action],
        'getargs': {'samtools_dir': ('download_samtools', 'dir')},
        'targets': [os.path.join(INSTALL_ROOT, 'bin', 'samtools')],
        'uptodate': [True]
    }


def task_install_bcftools():
    def action(bcftools_dir):
        sh('''
                make
                make prefix={} install
        '''.format(INSTALL_ROOT), cwd=bcftools_dir)

    return {
        'actions': [action],
        'task_dep': ['download_bcftools', 'download_htslib', 'install_zlib'],
        'getargs': {'bcftools_dir': ('download_bcftools', 'dir')},
        'targets': [os.path.join(INSTALL_BIN, 'bcftools')],
        'uptodate': [True]
    }


def task_install_bedtools():
    def action(bedtools_dir):
        sh('''
        make
        make prefix={} install
        '''.format(INSTALL_ROOT), cwd=bedtools_dir)

    return {
        'actions': [action],
        'getargs': {'bedtools_dir': ('download_bedtools', 'dir')},
        'task_dep': ['download_bedtools', 'install_zlib'],
        'targets': [os.path.join(INSTALL_BIN, 'bedtools')],
        'uptodate': [True]
    }


def task_install_gatk():
    def action(gatk_dir):
        sh('''
            GATK_JAR=`readlink -f target/GenomeAnalysisTK.jar`\
            cp $GATK_JAR {}
        '''.format(JAVA_LIBS_ROOT), cwd=gatk_dir)

    return {
        'actions': [action],
        'task_dep': ['download_maven', 'download_gatk'],
        'getargs': {'gatk_dir': ('download_gatk', 'dir')},
        'targets': [os.path.join(JAVA_LIBS_ROOT, 'GenomeAnalysisTK.jar')],
        'uptodate': [True]
    }


def task_install_perl_libs():
    """
    Installs all cpan libs from the cpan directory into the perl_lib directory
    :return:
    """

    def action(cpan_mirror_dir):
        # Use the cpan directory we made in download_perl_libs as a cpan mirror and install from there
        sh('cpanm -l {perl_lib} --mirror file://{cpan} --installdeps .'.format(perl_lib=PERL_LIB_ROOT,
                                                                                         cpan=cpan_mirror_dir),
           cwd=ROOT)

    return {
        'targets': [os.path.join(PERL_LIB_ROOT, 'bin')],
        'task_dep': ['install_perl', 'download_perl_libs', 'download_cpanm'],
        'actions': [action],
        'getargs': {'cpan_mirror_dir': ('download_perl_libs', 'dir')},
        'uptodate': [True]
    }
