import os
import re

BASE = os.path.realpath(os.path.join(__file__, '..', '..', '..'))
BATCHES = os.path.join(BASE, 'batches')
DESIGNS = os.path.join(BASE, 'designs')
CLASSPATH = os.path.join(BASE, 'tools/java_libs')
CONFIG_GROOVY = os.path.join(BASE, "pipeline", "config.groovy")
CONFIG_GROOVY_UTIL = os.path.join(BASE, "pipeline", "scripts", "config_groovy_util.sh")

def batch_dir(batch_name):
    return os.path.join(BATCHES, batch_name)

def read_config_groovy():
    """
    Parses the config groovy file and returns the interpolated values as a python dictionary
    """
    with open(CONFIG_GROOVY) as config_file:

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
        for key, val in result.iteritems():
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
