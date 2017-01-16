from pathlib import Path
import os

BASE = Path(os.environ.get('CPIPE_ROOT'))
BATCHES = BASE / 'batches'
DESIGNS = BASE / 'designs'
DESIGN_LIST = set([f.stem for f in DESIGNS.iterdir()])
CLASSPATH = BASE / 'tools/java_libs'
CONFIG_GROOVY = BASE / "pipeline" / "config.groovy"
CONFIG_GROOVY_UTIL = BASE / "pipeline" / "scripts" / "config_groovy_util.sh"
VERSION = BASE / 'version.txt'
PIPELINE_ID = BASE / 'pipeline_id'

