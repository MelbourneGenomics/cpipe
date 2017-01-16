import re
import pandas as pd
from typing import Union, Any

from .batch import Batch
from .design import Design
from .metadata import Metadata
from .paths import CONFIG_GROOVY, VERSION

def read_metadata(metadata_file: Any, parse_num=True):
    dtype = None if parse_num else str
    return pd.read_csv(metadata_file, sep='\t', dtype=dtype, na_values=[], keep_default_na=False)

def get_version():
    with open(VERSION) as version_file:
        return version_file.read().strip()

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


__all__ = ["Batch", "Design", "read_config_groovy", "read_metadata"]
