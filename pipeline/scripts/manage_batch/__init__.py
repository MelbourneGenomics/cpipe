#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path
from typing import List
from itertools import groupby

from cpipe_util.paths import CONFIG_GROOVY_UTIL, CLASSPATH, BASE, BATCHES, DESIGNS
from cpipe_util import read_metadata, Batch, Design, pathlib_patches
import cpipe_util
from manage_batch.schema import get_schema

def list_batches():
    df = cpipe_util.list_batches()
    df.to_csv(sys.stdout, sep='\t', index=False)



