"""
The root file for the doit build tool: http://pydoit.org/. This specifies tasks with dependencies and actions in the
 same way as a tool like make. Specifically, this file downloads and installs assets and build tools and runs other
 setup operations needed by cpipe. This module imports tasks from the tasks directory and executes whichever the user
 specified. The default task is 'install'.
"""
import sys

from tasks.common import has_swift_auth
from jinja2 import Template
from tasks.download import *
from tasks.install import *
from os import path
import re
import subprocess
from doit.tools import PythonInteractiveAction

from cpipe import get_version
from cpipe.scripts import create_bpipe_config
from multiprocessing import Lock

lock = Lock()

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
        'task_dep': ['copy_config', 'assets', 'generate_pipeline_id'],
        'clean': ['rm -rf {data} {tools} {tmp}/*'.format(data=DATA_ROOT, tools=TOOLS_ROOT, tmp=TMPDATA)]
    }


def task_assets():
    return {
        'actions': None,
        'task_dep': ['tool_assets', 'data_assets']
    }


def task_data_assets():
    return {
        'actions': None,
        'task_dep': [
            'download_dbnsfp',
            'install_vep_cache',
            'obtain_ucsc',
            'download_mills_and_1000g',
            'download_dbsnp',
            'obtain_trio_refinement',
            'download_chromosome_sizes',
            'download_vcfanno_data',
        ]
    }


def task_tool_assets():
    return {
        'actions': None,
        'task_dep': [
            'install_perl',
            'install_r',
            'install_bwa',
            'install_htslib',
            'install_samtools',
            'install_bcftools',
            'install_bedtools',
            'install_gatk',
            'install_vep',
            'install_fastqc',
            'install_groovy',
            'install_bpipe',
            'install_picard',
            'install_vep_libs',
            'install_vep_plugins',
            'install_java_libs',
            'install_bzip2',
            'install_xz',
            'install_pcre',
            'install_libcurl',
            'install_zlib',
            'install_vcfanno'
        ]
    }


def task_copy_config():
    return {
        'actions': None,
        'task_dep': ['copy_main_config'],
        'uptodate': [True]
    }


def task_copy_main_config():
    input = path.join(ROOT, 'pipeline', 'config.groovy.template')
    output = path.join(ROOT, 'pipeline', 'config.groovy')

    def action():
        with open(input, 'r') as input_file, open(output, 'w') as output_file:
            for line in input_file:
                substituted = line.replace("<ROOT DIR>", str(ROOT))
                output_file.write(substituted)

    return {
        'actions': [action],
        'targets': [output],
        'uptodate': [True]
    }


def task_generate_pipeline_id():
    """
    Creates a pipeline_id file, with an ID in the format HOSTNAME_2.X.X(version)_1
    :return:
    """
    def action(targets):
        version = get_version()

        for file in targets:
            with open(file, 'w') as id_file:
                id_file.write(f'{PIPELINE_ID}_{version}_1')

    return {
        'targets': ['pipeline_id'],
        'uptodate': [True],
        'actions': [action]
    }
