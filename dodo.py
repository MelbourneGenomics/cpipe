"""
The root file for the doit build tool: http://pydoit.org/. This specifies tasks with dependencies and actions in the
 same way as a tool like make. Specifically, this file downloads and installs assets and build tools and runs other
 setup operations needed by cpipe. This module imports tasks from the tasks directory and executes whichever the user
 specified. The default task is 'install'.
"""

from tasks.compile_tools import *
from tasks.common import ROOT, has_swift_auth
if has_swift_auth():
    from tasks.nectar_assets import *
else:
    from tasks.manual_assets import *
from os import path
import re
import subprocess

DOIT_CONFIG = {
    'default_tasks': ['install'],
    'backend': 'json'
}


def task_install():
    """
    Install cpipe in a manual install environment, e.g. not inside a docker container
    :return:
    """
    return {
        'actions': None,
        'task_dep': ['copy_config', 'assets', 'check_java']
    }


def task_assets():
    # If the user has all the necessary environment variables set, let them download the cached assets
    # from the object store
    if has_swift_auth():
        task = 'nectar_assets'
    else:
        task = 'manual_assets'

    return {
        'actions': None,
        'task_dep': [task]
    }


def task_copy_config():
    input = path.join(ROOT, 'pipeline', 'config.groovy.docker')
    output = path.join(ROOT, 'pipeline', 'config.groovy')

    def action():
        with open(input, 'r') as input_file, open(output, 'w') as output_file:
            for line in input_file:
                substituted = line.replace("<ROOT DIR>", ROOT)
                output_file.write(substituted)

    return {
        'actions': [action],
        'targets': [output],
        'uptodate': [True]
    }


def task_check_java():
    """
    Checks if we have a valid version of java. If we're in docker, automatically installs it, otherwise, throws
    an exception
    """

    def has_java_18():
        try:
            output = subprocess.check_output(
                "java -version",
                stderr=subprocess.STDOUT,
                shell=True
            )

            # Parse the major and minor versions to see if they have at least 1.8
            match = re.search(r'version "(?P<maj>\d+)\.(?P<min>\d+)\.(?P<build>[\d_]+)"', output)
            if int(match.groupdict()['maj']) >= 1 and int(match.groupdict()['min']) >= 8:
                return True
            else:
                return False
        except:
            return False

    def check_java():
        if has_java_18():
            return True
        else:
            raise OSError('Missing Java 1.8 or greater! Please install it to continue')

    return {
        'actions': [check_java],
        'task_dep': ['java_docker'] if in_docker() else []
    }
