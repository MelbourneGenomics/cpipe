import os
from pathlib import Path
import re
import pandas as pd
import sys

BASE = Path(__file__).parent.parent.parent
BATCHES = BASE / 'batches'
DESIGNS = BASE / 'designs'
DESIGN_LIST = set([f.stem for f in DESIGNS.iterdir()])
CLASSPATH = BASE / 'tools/java_libs'
CONFIG_GROOVY = BASE / "pipeline" / "config.groovy"
CONFIG_GROOVY_UTIL = BASE / "pipeline" / "scripts" / "config_groovy_util.sh"

def batch_dir(batch_name):
    return os.path.join(BATCHES, batch_name)

def list_batches(out):
    '''
        Prints the name of all batches that contain a samples.txt file
    '''

    # Find all directories that contain a samples.txt and add them to a list
    df = pd.DataFrame(columns=('Batch Name', 'Batch Path'))
    for root, dirs, files in os.walk(BATCHES):
        if 'samples.txt' in files:
            batch_name = os.path.basename(root)
            full_path = os.path.abspath(root)
            df = df.append({'Batch Name': batch_name, 'Batch Path': full_path}, ignore_index=True)

    # Sort them alphabetically by their batch name
    df = df.sort_values(by='Batch Name')

    # Print to output stream if requested, otherwise return the data frame
    if out is None:
        return df
    else:
        df.to_csv(out, sep='\t', index=False)

def read_config_groovy():
    """
    Parses the config groovy file and returns the interpolated values as a python dictionary
    """
    with CONFIG_GROOVY.open() as config_file:

        result = {}
        regex = re.compile(r'(?P<var>[\w_]+)=(?P<value>.+)')

        # Parse the config file one line at a time
        for line in config_file:

            # Ignore lines that are comments or blank
            if line.lstrip().startswith("//") or len(line) == 1:
                continue
            else:
                # Parse it using a regex
                match = regex.search(line).groupdict()
                key = match['var']
                val = match['value']

                # Give everything the correct data type
                if val == 'true':
                    val = True
                elif val == 'false':
                    val = False
                elif val.startswith('"'):
                    val = val.replace('"', '')
                else:
                    try:
                        val = float(val)
                    except:
                        pass

                result[key] = val

        # Now interpolate variables
        for key, val in list(result.items()):
            interpolate = re.compile(r'\$[\w_]+')
            if type(val) is str:
                while True:
                    match = interpolate.search(val)
                    if match:
                        var_interpolation = match.group(0)
                        var_name = var_interpolation[1:]
                        val = val.replace(var_interpolation, result[var_name])
                        result[key] = val
                    else:
                        break

        return result
