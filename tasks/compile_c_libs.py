"""
Tasks for compiling C libs that are needed by various runtimes (R, perl etc.)
All of these will try to install to the C_LIBS directory (using PREFIX) so that C_LIBS contains an
include, lib, bin etc. These will then be set to the path
"""

from tasks.common import *

def task_install_bzip2():
    def action(bzip2_dir):
        sh('''
            make -f Makefile-libbz2_so'
            make
            make install PREFIX={}
        '''.format(INSTALL_ROOT), cwd=bzip2_dir)
    return {
        'actions': [action],
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_bzip2'],
        'targets': [os.path.join(INSTALL_BIN, 'bzip2')],
        'uptodate': [True]
    }

def task_install_xz():
    def action(xz_dir):
        sh('''
            ./configure --prefix={}
            make
            make install
        '''.format(INSTALL_ROOT), cwd=xz_dir)
    return {
        'actions': [action],
        'task_dep': ['download_xz'],
        'targets': [os.path.join(INSTALL_BIN, 'xz')],
        'uptodate': [True]
    }

def task_install_pcre():
    def action(pcre_dir):
       sh('''
            ./configure --enable-utf8 --prefix={}
            make
            make install
       '''.format(INSTALL_ROOT), cwd=pcre_dir)
    return {
        'actions': [action],
        'task_dep': ['download_pcre'],
        'targets': [os.path.join(INSTALL_BIN, 'pcregrep')],
        'uptodate': [True]
    }



def task_install_libcurl():
    def action(libcurl_dir):
        sh('''
            ./configure --prefix={}
            make
            make install
        '''.format(INSTALL_ROOT), cwd=libcurl_dir)
    return {
        'actions': [action],
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_libcurl'],
        'targets': [os.path.join(INSTALL_BIN, 'curl')],
        'uptodate': [True]
    }


def task_install_zlib():
    def action(zlib_dir):
        sh('''
            ./configure --prefix={}
            make
            make install
        '''.format(INSTALL_ROOT), cwd=zlib_dir)
    return {
        'actions': [action],
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_zlib'],
        'targets': [os.path.join(INSTALL_BIN, 'libz.so')],
        'uptodate': [True]
    }


