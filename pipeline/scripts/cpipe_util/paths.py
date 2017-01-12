from pathlib import Path


BASE = Path(__file__).parent.parent.parent.parent.resolve()
BATCHES = BASE / 'batches'
DESIGNS = BASE / 'designs'
DESIGN_LIST = set([f.stem for f in DESIGNS.iterdir()])
CLASSPATH = BASE / 'tools/java_libs'
CONFIG_GROOVY = BASE / "pipeline" / "config.groovy"
CONFIG_GROOVY_UTIL = BASE / "pipeline" / "scripts" / "config_groovy_util.sh"

