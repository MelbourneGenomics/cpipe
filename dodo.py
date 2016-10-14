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
        'task_dep': ['copy_config', 'assets'],
        'clean': ['rm -rf {data} {tools} {tmp}/*'.format(data=DATA_ROOT, tools=TOOLS_ROOT, tmp=TMPDATA)]
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
    input = path.join(ROOT, 'pipeline', 'config.groovy.template')
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

    def parse_java_version():
        """
        Parses the Java version and returns a tuple of (maj,min,build)
        :return:
        """
        output = subprocess.check_output(
            "java -version",
            stderr=subprocess.STDOUT,
            shell=True
        )

        # Parse the major and minor versions using a regex
        parsed = re.search(r'version "(?P<maj>\d+)\.(?P<min>\d+)\.(?P<build>[\d_]+)"', output).groupdict()
        return (int(parsed['maj']), int(parsed['min']), parsed['build'])

    def check_java():
        (maj, min, build) = parse_java_version()

        if not maj >= 1 and min >= 8:
            raise OSError('Missing Java 1.8 or greater! Please install it to continue')

        if build == '0_20':
            raise OSError(
                'You have Java 1.8 Build 20! This version of Java has compatibility issues with groovy bytecode.'
                ' Please either upgrade your version of java, or, if this is not possible, edit your pipeline/config.groovy'
                ' file and replace \'JAVA="java"\' with \'JAVA="java -noverify"\', then re-run this script with'
                ' the --no-java-check flag.'
            )

    return {
        'actions': [check_java],
        'task_dep': ['java_docker'] if in_docker() else []
    }
