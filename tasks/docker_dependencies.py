"""
The purpose of tasks in here is that various tasks, e.g. `compile_r` should check if they're in docker, and if they are,
add one of these tasks e.g. `task_r_docker_dependencies` as a task_dep. If they're not in docker, the compile step will
go ahead, and hopefully succeed, but if it fails (e.g. missing a dependency), it will alert the user and they will have
to install the relevant headers themselves
"""

def task_r_docker_dependencies():
    return {
        'actions': ['apt-get install -y gfortran xorg-dev libreadline-dev libbz2-dev liblzma-dev libpcre3-dev'],
        'uptodate': [False]
    }
