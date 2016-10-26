import os
import glob
from tasks.nectar.nectar_util import *

from tasks.install.install_c_libs import *


def task_install_perl():
    def action(perl_dir):
        sh('''
              ./Configure -de -Dprefix={}
              make
              make install
        '''.format(INSTALL_ROOT),
           cwd=perl_dir
           )
        if has_swift_auth():
            add_to_manifest('perl')

    return {
        'actions': [action],
        'getargs': {'perl_dir': ('download_perl', 'dir')},
        'task_dep': ['download_perl'],
        'targets': [os.path.join(INSTALL_ROOT, 'bin', 'perl')],
        'uptodate': [not nectar_asset_needs_update('perl')],
    }


def task_install_r():
    def action(r_dir):
        sh('''
            ./configure
             make
             make prefix={} install
        '''.format(r_dir), cwd=r_dir)
        if has_swift_auth():
            add_to_manifest('r')

    return {
        'actions': [action],
        'task_dep': ['download_perl', 'download_r', 'install_bzip2', 'install_xz', 'install_pcre', 'install_libcurl',
                     'install_zlib'],
        'getargs': {'r_dir': ('download_r', 'dir')},
        'targets': [os.path.join(INSTALL_ROOT, 'bin', 'R')],
        'uptodate': [not nectar_asset_needs_update('r')],
    }


def task_install_bwa():
    def action(bwa_dir):
        sh('make', cwd=bwa_dir)
        if has_swift_auth():
            add_to_manifest('bwa')

    return {
        'actions': [action],
        'getargs': {'bwa_dir': ('download_bwa', 'dir')},
        'task_dep': ['download_bwa'],
        'targets': [os.path.join(INSTALL_BIN, 'bwa')],
        'uptodate': [not nectar_asset_needs_update('bwa')],
    }


def task_install_htslib():
    def action(htslib_dir):
        sh('''
        ./configure --prefix={0}
        make
        make prefix={0} install
        '''.format(INSTALL_ROOT), cwd=htslib_dir)
        if has_swift_auth():
            add_to_manifest('htslib')

    return {
        'task_dep': ['download_htslib', 'install_zlib'],
        'actions': [action],
        'getargs': {'htslib_dir': ('download_htslib', 'dir')},
        'targets': [os.path.join(INSTALL_ROOT, 'bin', 'htsfile')],
        'uptodate': [not nectar_asset_needs_update('htslib')],
    }


def task_install_samtools():
    def action(samtools_dir):
        sh('''
                ./configure --prefix={0} --with-htslib={1}
                make
                make prefix={0} install
        '''.format(INSTALL_ROOT, INSTALL_LIB), cwd=samtools_dir)
        if has_swift_auth():
            add_to_manifest('samtools')

    return {
        'task_dep': ['install_zlib', 'install_htslib', 'download_samtools'],
        'actions': [action],
        'getargs': {'samtools_dir': ('download_samtools', 'dir')},
        'targets': [os.path.join(INSTALL_ROOT, 'bin', 'samtools')],
        'uptodate': [not nectar_asset_needs_update('samtools')],
    }


def task_install_bcftools():
    def action(bcftools_dir):
        sh('''
                make
                make prefix={} install
        '''.format(INSTALL_ROOT), cwd=bcftools_dir)
        if has_swift_auth():
            add_to_manifest('bcftools')

    return {
        'actions': [action],
        'task_dep': ['download_bcftools', 'download_htslib', 'install_zlib'],
        'getargs': {'bcftools_dir': ('download_bcftools', 'dir')},
        'targets': [os.path.join(INSTALL_BIN, 'bcftools')],
        'uptodate': [not nectar_asset_needs_update('bcftools')],
    }


def task_install_bedtools():
    def action(bedtools_dir):
        sh('''
        make
        make prefix={} install
        '''.format(INSTALL_ROOT), cwd=bedtools_dir)
        if has_swift_auth():
            add_to_manifest('bedtools')

    return {
        'actions': [action],
        'getargs': {'bedtools_dir': ('download_bedtools', 'dir')},
        'task_dep': ['download_bedtools', 'install_zlib'],
        'targets': [os.path.join(INSTALL_BIN, 'bedtools')],
        'uptodate': [not nectar_asset_needs_update('bedtools')],
    }


def task_install_gatk():
    def action(gatk_dir):
        sh('''
            GATK_JAR=`readlink -f target/GenomeAnalysisTK.jar`\
            cp GenomeAnalysisTK.jar {}
        '''.format(JAVA_LIBS_ROOT), cwd=gatk_dir)
        if has_swift_auth():
            add_to_manifest('gatk')

    return {
        'actions': [action],
        'task_dep': ['download_maven', 'download_gatk'],
        'getargs': {'gatk_dir': ('download_gatk', 'dir')},
        'targets': [os.path.join(JAVA_LIBS_ROOT, 'GenomeAnalysisTK.jar')],
        'uptodate': [not nectar_asset_needs_update('gatk')],
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
        if has_swift_auth():
            add_to_manifest('perl_libs')

    return {
        'targets': [os.path.join(PERL_LIB_ROOT, 'bin')],
        'task_dep': ['install_perl', 'download_perl_libs', 'download_cpanm'],
        'actions': [action],
        'getargs': {'cpan_mirror_dir': ('download_perl_libs', 'dir')},
        'uptodate': [not nectar_asset_needs_update('perl_libs')],
    }


def task_install_vep():
    def action(vep_dir):
        delete_and_copy(vep_dir, VEP_ROOT)
        if has_swift_auth():
            add_to_manifest('vep')
    return {
        'actions': [action],
        'targets': [VEP_ROOT, os.path.join(VEP_ROOT, 'variant_effect_predictor.pl')],
        'uptodate': [not nectar_asset_needs_update('vep')],
        'getargs': {'vep_dir': ('download_vep', 'dir')},
    }

def task_install_fastqc():
    script_bin = os.path.join(INSTALL_BIN, 'fastqc')

    def action(fastqc_dir):
        delete_and_copy(fastqc_dir, FASTQC_ROOT)
        os.symlink(os.path.join(FASTQC_ROOT, 'fastqc'), script_bin)
        if has_swift_auth():
            add_to_manifest('fastqc')

    return {
        'actions': [action],
        'targets': [script_bin, FASTQC_ROOT],
        'uptodate': [not nectar_asset_needs_update('fastqc')],
        'getargs': {'fastqc_dir': ('download_fastqc', 'dir')},
    }

def task_install_bpipe():
    def action(bpipe_dir):
        delete_and_copy(bpipe_dir, BPIPE_ROOT)
        if has_swift_auth():
            add_to_manifest('bpipe')

    return {
        'actions': [action],
        'targets': [BPIPE_ROOT, os.path.join(BPIPE_ROOT, 'bin', 'bpipe')],
        'uptodate': [not nectar_asset_needs_update('bpipe')],
        'getargs': {'bpipe_dir': ('download_bpipe', 'dir')},
    }

def task_install_picard():
    picard_target = os.path.join(JAVA_LIBS_ROOT, 'picard.jar')
    def action(picard_dir):
        picard_jar = os.path.join(picard_dir, 'picard.jar')
        delete_and_copy(picard_jar, picard_target)
        if has_swift_auth():
            add_to_manifest('picard')

    return {
        'actions': [action],
        'targets': [picard_target],
        'uptodate': [not nectar_asset_needs_update('picard')],
        'getargs': {'picard_dir': ('download_picard', 'dir')},
    }

def task_install_vep_libs():
    def action(vep_libs_dir):
        delete_and_copy(vep_libs_dir, VEP_LIBS_ROOT)
        if has_swift_auth():
            add_to_manifest('vep_libs')

    return {
        'actions': [action],
        'uptodate': [not nectar_asset_needs_update('vep_libs')],
        'targets': [VEP_LIBS_ROOT, os.path.join(VEP_LIBS_ROOT, 'Bio', 'TreeIO.pm')],
        'getargs': {'vep_libs_dir': ('download_vep_libs', 'dir')},
    }

def task_install_vep_plugins():
    def action(vep_plugins_dir):
        delete_and_copy(vep_plugins_dir, VEP_PLUGIN_ROOT)
        if has_swift_auth():
            add_to_manifest('vep_plugins')

    return {
        'actions': [action],
        'uptodate': [not nectar_asset_needs_update('vep_plugins')],
        'targets': [VEP_PLUGIN_ROOT, os.path.join(VEP_PLUGIN_ROOT, 'GO.pm')],
        'getargs': {'vep_plugins_dir': ('download_vep_plugins', 'dir')},
    }

def task_install_junit_xml_formatter():
    target = os.path.join(JAVA_LIBS_ROOT, 'JUnitXmlFormatter.jar')
    def action(junit_xml_dir):
        jar = glob.glob(os.path.join(junit_xml_dir, 'JUnitXmlFormatter*'))[0]
        delete_and_copy(jar, target)
        if has_swift_auth():
            add_to_manifest('junit_xml_formatter')

    return {
        'actions': [action],
        'targets': [target],
        'uptodate': [not nectar_asset_needs_update('junit_xml_formatter')],
        'getargs': {'junit_xml_dir': ('download_junit_xml_formatter', 'dir')},
    }

def task_install_groovy_ngs_utils():
    target = os.path.join(JAVA_LIBS_ROOT, 'groovy-ngs-utils.jar')
    def action(groovy_ngs_dir):
        jar = os.path.join(groovy_ngs_dir, 'build' 'libs' 'groovy-ngs-utils.jar')
        delete_and_copy(jar, target)
        if has_swift_auth():
            add_to_manifest('groovy_ngs_utils')

    return {
        'actions': [action],
        'targets': [target],
        'uptodate': [not nectar_asset_needs_update('groovy_ngs_utils')],
        'getargs': {'groovy_ngs_dir': ('download_groovy_ngs_utils', 'dir')},
    }

def task_install_takari_cpsuite():
    target = os.path.join(JAVA_LIBS_ROOT, 'takari-cpsuite.jar')
    def action(takari_cpsuite_dir):
        jar = glob.glob(os.path.join(takari_cpsuite_dir, 'takari-cpsuite*.jar'))[0]
        delete_and_copy(jar, target)
        if has_swift_auth():
            add_to_manifest('takari_cpsuite')

    return {
        'actions': [action],
        'targets': [target],
        'uptodate': [not nectar_asset_needs_update('takari_cpsuite')],
        'getargs': {'takari_cpsuite_dir': ('download_takari_cpsuite', 'dir')},
    }

def task_install_maven():
    target = os.path.join(MAVEN_ROOT, 'bin', 'mvn')
    def action(maven_dir):
        delete_and_copy(maven_dir, MAVEN_ROOT)
        if has_swift_auth():
            add_to_manifest('maven')

    return {
        'actions': [action],
        'targets': [target],
        'uptodate': [True],
        'getargs': {'maven_dir': ('download_maven', 'dir')},
    }
