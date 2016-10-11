from tasks.common import *

def task_compile_bzip2():
    return {
        'actions': [
            cmd('make -f Makefile-libbz2_so', cwd=BZIP_ROOT),
            cmd('make', cwd=BZIP_ROOT),
        ],
        'task_dep': ['download_nectar_assets' if has_swift_auth() else 'download_bzip2'],
        'targets': [os.path.join(BZIP_ROOT, 'bzip2')],
        'uptodate': [True]
    }

