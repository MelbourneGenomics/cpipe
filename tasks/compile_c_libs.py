"""
Tasks for compiling C libs that are needed by various runtimes (R, perl etc.)
All of these will try to install to the C_LIBS directory (using PREFIX) so that C_LIBS contains an
include, lib, bin etc. These will then be set to the path
"""

from tasks.common import *

def task_compile_bzip2():
    return {
        'actions': [
            cmd('make -f Makefile-libbz2_so', cwd=BZIP_ROOT),
            cmd('make', cwd=BZIP_ROOT),
            cmd('make install PREFIX={}'.format(C_INCLUDE_ROOT), cwd=BZIP_ROOT),
        ],
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_bzip2'],
        'targets': [os.path.join(C_INCLUDE_ROOT, 'bin', 'bzip2')],
        'uptodate': [True]
    }

def task_compile_xz():
    return {
        'actions': [
            cmd('./configure --prefix={}'.format(C_INCLUDE_ROOT), cwd=XZ_ROOT),
            cmd('make', cwd=XZ_ROOT),
            cmd('make install'.format(C_INCLUDE_ROOT), cwd=XZ_ROOT),
        ],
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_xz'],
        'targets': [os.path.join(C_INCLUDE_ROOT, 'bin', 'xz')],
        'uptodate': [True]
    }

def task_compile_pcre():
    return {
        'actions': [
            cmd('./configure --enable-utf8 --prefix={}'.format(C_INCLUDE_ROOT), cwd=PCRE_ROOT),
            cmd('make', cwd=PCRE_ROOT),
            cmd('make install'.format(C_INCLUDE_ROOT), cwd=PCRE_ROOT),
        ],
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_pcre'],
        'targets': [os.path.join(C_INCLUDE_ROOT, 'bin', 'pcregrep')],
        'uptodate': [True]
    }


