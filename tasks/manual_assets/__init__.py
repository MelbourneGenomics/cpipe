"""
Defines the doit tasks required to download cpipe assets without access to NECTAR. Must be run through the root doit.py
script
"""

from tasks.manual_assets.dependencies import *
from tasks.manual_assets.download_reference_files import *
from tasks.manual_assets.download_tools import *

def task_manual_assets():
    return {
        'actions': None,
        'task_dep': ['tool_assets', 'data_assets', 'compile_all']
    }