import subprocess
import re

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
    parsed = re.search(br'version "(?P<maj>\d+)\.(?P<min>\d+)\.(?P<build>[\d_]+)"', output).groupdict()
    return int(parsed['maj']), int(parsed['min']), parsed['build']

def check_java():
    """
    Checks if we have a valid version of java. If we're in docker, automatically installs it, otherwise, throws
    an exception
    """
    (maj, min, build) = parse_java_version()

    if not (maj >= 1 and min >= 8):
        raise OSError('Missing Java 1.8 or greater! Please install it to continue')

    if maj == 1 and min == 8 and build == '0_20':
        raise OSError(
            'You have Java 1.8 build 20! This version of Java has compatibility issues with groovy bytecode.'
            ' Please either upgrade your version of java, or, if this is not possible, edit your pipeline/config.groovy'
            ' file and set the JAVA_OPTS variable to \'JAVA_OPTS="-noverify"\', then re-run this script with'
            ' the --no-java-check flag.'
        )

def main():
    check_java()

if __name__ == '__main__':
    main()

