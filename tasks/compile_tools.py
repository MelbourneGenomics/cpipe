from tasks.common import *
from tasks.docker_dependencies import *
from doit.tools import create_folder
from tasks.compile_c_libs import *
import os
import subprocess

def task_compile_all():
    'Compiles all assets including java assets; this is only needed for a manual install'
    return {
        'actions': None,
        'task_dep': [
            'compile_perl',
            'compile_r',
            'compile_bwa',
            'compile_htslib',
            'compile_samtools',
            'compile_bcftools',
            'compile_bedtools',
            'compile_gatk',
            'install_perl_libs'
        ]
    }

def task_compile_nectar():
    'Compiles assets downloaded from nectar; excludes java compilation since this should already have been done'
    return {
        'actions': None,
        'task_dep': [
            'compile_perl',
            'compile_r',
            'compile_bwa',
            'compile_htslib',
            'compile_samtools',
            'compile_bcftools',
            'compile_bedtools', 
            'install_perl_libs'
        ]
    }

def task_compile_perl():

    return {
        'actions': [
            cmd('yes | ./configure.sh && make', cwd=PERL_ROOT)
        ],
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_perl'],
        'targets': [os.path.join(PERL_ROOT, 'perl')],
        'uptodate': [True]
    }

def task_compile_r():
    task_dep = ['download_nectar_assets'] if has_swift_auth() else ['download_perl', 'download_r']
    task_dep.append('compile_bzip2')
    task_dep.append('compile_xz')
    task_dep.append('compile_pcre')
    task_dep.append('compile_libcurl')
    task_dep.append('compile_zlib')
    if in_docker():
        task_dep.append('r_docker_dependencies')

    def action(): 
        subprocess.check_call('source {env} && ./configure && make'.format(env=ENVIRONMENT_FILE), cwd=R_ROOT, env=get_c_env(), shell=True, executable='bash')

    return {
        'actions': [action],
        'task_dep': task_dep,
        'targets': [os.path.join(R_ROOT, 'bin', 'R')],
        'uptodate': [True]
    }

def task_compile_bwa():
    return {
        'actions': [
            cmd('make', cwd=BWA_ROOT)
        ],
        'task_dep': [
            'download_nectar_assets' if has_swift_auth() else 'download_bwa'],
        'targets': [os.path.join(BWA_ROOT, 'bwa')],
        'uptodate': [True]
    }

def task_compile_htslib():
    return {
        'actions': [
            cmd('make', cwd=HTSLIB_ROOT)
        ],
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_htslib'],
        'targets': [os.path.join(HTSLIB_ROOT, 'htsfile')],
        'uptodate': [True]
    }

def task_compile_samtools():
    task_dep = ['download_nectar_assets'] if has_swift_auth() else ['download_samtools', 'download_htslib']

    if in_docker():
        task_dep.append('samtools_docker_dependencies')

    return {
        'actions': [
            cmd('make', cwd=SAMTOOLS_ROOT)
        ],
        'task_dep': task_dep,
        'targets': [os.path.join(SAMTOOLS_ROOT, 'samtools')],
        'uptodate': [True]
    }

def task_compile_bcftools():
    return {
        'actions': [
            cmd('make', cwd=BCFTOOLS_ROOT)
        ],
        'task_dep': ['download_nectar_assets'] if has_swift_auth() else ['download_bcftools', 'download_htslib'],
        'targets': [os.path.join(BCFTOOLS_ROOT, 'bcftools')],
        'uptodate': [True]
    }

def task_compile_bedtools():
    task_dep = ['download_nectar_assets' if has_swift_auth() else 'download_bedtools']
    task_dep.append('compile_zlib')

    return {
        'actions': [
            cmd('make', cwd=BEDTOOLS_ROOT)
        ],
        'task_dep': task_dep,
        'targets': [os.path.join(BEDTOOLS_ROOT, 'bin', 'bedtools')],
        'uptodate': [True]
    }

def task_compile_gatk():
    return {
        'actions': [
            cmd('''
            mvn verify -P\!queue\
            && GATK_JAR=`readlink -f target/GenomeAnalysisTK.jar`\
            && unlink target/GenomeAnalysisTK.jar\
            && mv $GATK_JAR ./GenomeAnalysisTK.jar\
            && bash -O extglob -c 'rm -rf !(GenomeAnalysisTK.jar)'
            ''', cwd=GATK_ROOT)
        ],
        'task_dep': ['download_nectar_assets'] if has_swift_auth() else ['download_maven', 'download_gatk'],
        'targets': [os.path.join(GATK_ROOT, 'GenomeAnalysisTK.jar')],
        'uptodate': [True]
    }


def task_install_perl_libs():
    """
    Installs all cpan libs from the cpan directory into the perl_lib directory
    :return:
    """

    return {
        'targets': [os.path.join(PERL_LIB_ROOT, 'bin')],
        'task_dep': ['download_nectar_assets'] if has_swift_auth() else ['download_perl_libs', 'download_cpanm'],
        'actions': [
            # Use the cpan directory we made in download_perl_libs as a cpan mirror and install from there
            cmd('cpanm -l {perl_lib} --mirror file://{tools_dir}/cpan --installdeps .'.format(tools_dir=TOOLS_ROOT, perl_lib=PERL_LIB_ROOT), cwd=ROOT, env=get_cpanm_env())
        ],
        'uptodate': [True]
    }
