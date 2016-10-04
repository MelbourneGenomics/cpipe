from tasks.common import *
from tasks.docker_dependencies import *
from doit.tools import create_folder

def task_compile_tools():
    return {
        'actions': [None],
        'task_dep': [
            'compile_perl',
            'compile_r',
            'compile_bwa',
            'compile_htslib',
            'compile_samtools',
            'compile_bcftools',
            'compile_bedtools',
            'compile_gatk'
        ]
    }

def task_compile_perl():
    return {
        'actions': [
            cmd('yes | ./configure.sh && make', cwd=PERL_ROOT)
        ],
        'task_dep': ['download_perl'],
        'targets': [os.path.join(PERL_ROOT, 'perl')],
        'uptodate': [True]
    }

def task_compile_r():
    task_dep = ['download_r']
    if in_docker():
        task_dep.append('r_docker_dependencies')

    return {
        'actions': [
            cmd('./configure && make', cwd=R_ROOT)
        ],
        'task_dep': task_dep,
        'targets': [os.path.join(R_ROOT, 'bin', 'R')],
        'uptodate': [True]
    }

def task_compile_bwa():
    return {
        'actions': [
            cmd('./configure && make', cwd=BWA_ROOT)
        ],
        'task_dep': ['download_bwa'],
        'targets': [os.path.join(BWA_ROOT, 'bwa')],
        'uptodate': [True]
    }

def task_compile_htslib():
    return {
        'actions': [
            cmd('make', cwd=HTSLIB_ROOT)
        ],
        'task_dep': ['download_htslib'],
        'targets': [os.path.join(HTSLIB_ROOT, 'htsfile')],
        'uptodate': [True]
    }

def task_compile_samtools():
    return {
        'actions': [
            cmd('make', cwd=SAMTOOLS_ROOT)
        ],
        'task_dep': ['download_samtools'],
        'targets': [os.path.join(SAMTOOLS_ROOT, 'samtools')],
        'uptodate': [True]
    }

def task_compile_bcftools():
    return {
        'actions': [
            cmd('make', cwd=BCFTOOLS_ROOT)
        ],
        'task_dep': ['download_bcftools'],
        'targets': [os.path.join(BCFTOOLS_ROOT, 'bcftools')],
        'uptodate': [True]
    }

def task_compile_bedtools():
    return {
        'actions': [
            cmd('make', cwd=BEDTOOLS_ROOT)
        ],
        'task_dep': ['download_bedtools'],
        'targets': [os.path.join(BEDTOOLS_ROOT, 'bedtools')],
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
        'task_dep': ['download_gatk'],
        'targets': [os.path.join(GATK_ROOT, 'gatk')],
        'uptodate': [True]
    }


def task_install_perl_libs():
    """
    Installs all cpan libs from the cpan directory into the perl_lib directory
    :return:
    """
    return {
        'targets': [PERL_LIB_ROOT],
        'task_dep': ['download_perl_libs'],
        'actions': [
            lambda: create_folder(PERL_LIB_ROOT),
            cmd('cpanm --mirror file://{0}/cpan -L {0}/perl_lib --installdeps .'.format(TOOLS_ROOT), cwd=INSTALL_ROOT)
        ],
        # 'uptodate': [True]
    }
