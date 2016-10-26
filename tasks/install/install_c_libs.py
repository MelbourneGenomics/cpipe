"""
Tasks for compiling C libs that are needed by various runtimes (R, perl etc.)
All of these will try to install to the C_LIBS directory (using PREFIX) so that C_LIBS contains an
include, lib, bin etc. These will then be set to the path
"""

from tasks.common import *
from tasks.nectar.nectar_util import *

def task_install_bzip2():
    def action(bzip2_dir):
        sh('''
            make -f Makefile-libbz2_so
            make
            make install PREFIX={}
        '''.format(INSTALL_ROOT), cwd=bzip2_dir)
        if has_swift_auth():
            add_to_manifest('bzip2')
    return {
        'actions': [action],
        'setup': ['download_bzip2'],
        'targets': [os.path.join(INSTALL_BIN, 'bzip2')],
        'uptodate': [not nectar_asset_needs_update('bzip2')],
        'getargs': {'bzip2_dir': ('download_bzip2', 'dir')},
    }

def task_install_xz():
    def action(xz_dir):
        sh('''
            ./configure --prefix={}
            make
            make install
        '''.format(INSTALL_ROOT), cwd=xz_dir)
        if has_swift_auth():
            add_to_manifest('xz')
    return {
        'actions': [action],
        'targets': [os.path.join(INSTALL_BIN, 'xz')],
        'setup': ['download_xz'],
        'uptodate': [not nectar_asset_needs_update('xz')],
        'getargs': {'xz_dir': ('download_xz', 'dir')},
    }

def task_install_pcre():
    def action(pcre_dir):
       sh('''
            ./configure --enable-utf8 --prefix={}
            make
            make install
       '''.format(INSTALL_ROOT), cwd=pcre_dir)
       if has_swift_auth():
           add_to_manifest('pcre')

    return {
        'actions': [action],
        'setup': ['download_pcre'],
        'targets': [os.path.join(INSTALL_BIN, 'pcregrep')],
        'uptodate': [not nectar_asset_needs_update('pcre')],
        'getargs': {'pcre_dir': ('download_pcre', 'dir')},
    }



def task_install_libcurl():
    def action(libcurl_dir):
        sh('''
            ./configure --prefix={}
            make
            make install
        '''.format(INSTALL_ROOT), cwd=libcurl_dir)
        if has_swift_auth():
            add_to_manifest('libcurl')
    return {
        'actions': [action],
        'setup': ['download_libcurl'],
        'targets': [os.path.join(INSTALL_BIN, 'curl')],
        'uptodate': [not nectar_asset_needs_update('libcurl')],
        'getargs': {'libcurl_dir': ('download_libcurl', 'dir')},
    }


def task_install_zlib():
    def action(zlib_dir):
        sh('''
            ./configure --prefix={}
            make
            make install
        '''.format(INSTALL_ROOT), cwd=zlib_dir)
        if has_swift_auth():
            add_to_manifest('zlib')
    return {
        'actions': [action],
        'targets': [os.path.join(INSTALL_ROOT, 'lib', 'libz.so')],
        'setup': ['download_zlib'],
        'uptodate': [not nectar_asset_needs_update('zlib')],
        'getargs': {'zlib_dir': ('download_zlib', 'dir')},
    }


