"""
The root file for the doit build tool: http://pydoit.org/. This specifies tasks with dependencies and actions in the
 same way as a tool like make. Specifically, this file downloads and installs assets and build tools and runs other
 setup operations needed by cpipe. This module imports tasks from the tasks directory and executes whichever the user
 specified. The default task is 'install'.
"""

from tasks.manual_assets import *
from tasks.nectar_assets import *
from tasks.common import ROOT

import os

DOIT_CONFIG = {'default_tasks': ['install']}

def task_install():
    """
    Install cpipe in a manual install environment, e.g. not inside a docker container
    :return:
    """
    return {
        'actions': [None],
        'task_dep': ['assets', 'copy_config']
    }

def task_docker_install():
    """
    Install cpipe for a docker container
    :return:
    """
    return {
        'actions': [None],
        'task_dep': ['assets', 'copy_docker_config']
    }

def task_assets():
    swift_credentials = {
        'OS_AUTH_URL',
        'OS_TENANT_ID'
        'OS_TENANT_NAME',
        'OS_PROJECT_NAME',
        'OS_USERNAME',
        'OS_PASSWORD'
    }

    # If the user has all the necessary environment variables set, let them download the cached assets
    # from the object store
    if swift_credentials.issubset(os.environ.keys()):
        task = 'nectar_assets'
    else:
        task = 'manual_assets'

    return {
        'actions': [None],
        'task_dep': [task]
    }

def task_copy_config():
    return {
        'actions': [lambda: shutil.copy(
            path.join(ROOT, 'pipeline', 'config.groovy.template'),
            path.join(ROOT, 'pipeline', 'config.groovy')
        )],
    }

def task_copy_docker_config():
    return {
        'actions': [lambda: shutil.copy(
            path.join(ROOT, 'pipeline', 'config.groovy.docker'),
            path.join(ROOT, 'pipeline', 'config.groovy')
        )]
    }
