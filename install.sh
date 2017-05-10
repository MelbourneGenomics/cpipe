#!/usr/bin/env bash

# Fail on error
set -e

# Set variables
PYTHON_VERSION='3.6.0'
PYTHON_INTERPRETER='python3.6'
ROOT=$(readlink -f $(dirname ${BASH_SOURCE}))
export TMPDIR=${ROOT}/tmpdata # Write temporary files to tmpdata
export CPIPE_ROOT=$ROOT
TEMP_SUBDIR=`mktemp -d`
SYS_PYTHON=${ROOT}/tools
SYS_PYBIN=${SYS_PYTHON}/bin
SYS_INTERPRETER=${SYS_PYBIN}/${PYTHON_INTERPRETER}
PYTHON=${ROOT}/tools/python
VENV=${PYTHON}/bin/activate

# Printing utilities
if [[ $- == *i* ]]; then
    bold=$(tput bold)
    normal=$(tput sgr0)
else
    bold=
    normal=
fi

# Usage function
function usage {
  echo "${bold}Cpipe Installer"
  echo "  ${bold}--help, --usage"
  echo "    ${normal}Print this help page to stdout"
  echo "  ${bold}-n, --processes <process number>"
  echo "    ${normal}Set the maximum number of processes to use for the install. The higher number the faster the install, but the more memory used. Defaults to the output of 'nproc --all', the number of available processing units (currently `nproc --all` on your system)"
  echo "  ${bold}-c, --credentials </path/to/swift_credentials.sh>"
  echo "    ${normal}Use the specified swift credentials file to download assets from NECTAR. Defaults to looking in the cpipe root directory"
  echo "  ${bold}-v, --verbose"
  echo "    ${normal}Print everything to stdout instead of just errors. Good for debugging an install"
  echo "  ${bold}-s, --no-swift"
  echo "    ${normal}Do a manual install instead of downloading assets from NECTAR. Strongly NOT recommended as this will potentially take days to complete"
  echo "  ${bold}-t, --task <taskname>"
  echo "    ${normal}Specify one or more tasks to run instead of a full install. Each time you use this flag it will add another task. Don't use this unless you know what you're doing"
  echo "  ${bold}-p, --no-pip"
  echo "    ${normal}Don't update pip modules. Don't use this unless you know what you're doing"
}

# Parse arguments
ARGS=$(getopt -o n:pvt:c:s --long "processes:,verbose,task:,no-pip,credentials:,help,usage,no-swift" -n $(basename $0) -- "$@")
eval set -- "$ARGS"

PROCESSES=`nproc --all`
VERBOSITY=1
USE_PIP=1
TASKS='install'
CUSTOM_TASKS=''
CREDENTIALS="${ROOT}/swift_credentials.sh"
USE_SWIFT=1
MODE='auto'

while true ; do
    case "$1" in
        -n|--processes)
          PROCESSES=$2
          shift 2;;
        -v|--verbose)
          VERBOSITY=2
          shift 1 ;;
        -p|--no-pip)
          USE_PIP=0
          shift 1 ;;
        --usage|--help)
          usage
          exit 0;;
        -t|--task)
          CUSTOM_TASKS="${CUSTOM_TASKS} $2"
          shift 2;;
        -c|--credentials)
          CREDENTIALS=$2
          shift 2;;
        -s|--no-swift)
          USE_SWIFT=0
          shift 1;;
        --)
          break ;;
        *)
          >&2 echo "Invalid argument \"$1\""
          usage
          exit 1 ;;
    esac
done

# If the user specified any tasks, do them instead of install
if [[ -n $CUSTOM_TASKS ]] ; then
    TASKS=$CUSTOM_TASKS
fi

# If credentials were specified, load them
if (( USE_SWIFT )); then

    # If the credentials exist, load them. Otherwise, warn the user they don't have credentials
    if [[ ! -f $CREDENTIALS ]]; then
        >&2 echo "Credentials file doesn't exist at the path specified ($CREDENTIALS). Cpipe will perform an incomplete\
 installation. If you are part of MGHA you probably want to obtain download credentials and place them in the path above."
    else
        >&2 echo "Credentials file found in path specified ($CREDENTIALS). Performing full install."
        # Load swift credentials
        source ${ROOT}/swift_credentials.sh
    fi
else
    MODE='manual'
fi

# Output stream
if [[ $VERBOSITY == 2 ]] ; then
    OUTPUT_STREAM="/dev/stdout"
else
    OUTPUT_STREAM="/dev/null"
fi

    # Use python-build to install python
    if [[ ! -f ${VENV} ]]; then

            echo -n 'Installing local python...'

            {
                pushd ${TEMP_SUBDIR}
                        git clone --depth 1 git://github.com/yyuu/pyenv.git
                        pyenv/plugins/python-build/bin/python-build ${PYTHON_VERSION} ${SYS_PYTHON}
                popd

                # Make a virtual environment
                ${SYS_INTERPRETER} -m venv ${PYTHON}

                # Delete the temporary files
                rm -rf ${TEMP_SUBDIR}

            } > ${OUTPUT_STREAM}

    fi

    {
        # Load virtualenv
        source ${ROOT}/_env

        # Install pip dependencies
        if (( USE_PIP )); then
            pip install --upgrade setuptools pip
            pip install -e ${ROOT}/lib -q
        fi ;

    } > ${OUTPUT_STREAM}

# Run the interactive task first
doit --verbosity $VERBOSITY copy_bpipe_config mode=${MODE}

# Now run the full install, which is all automated
doit -n $PROCESSES --verbosity $VERBOSITY $TASKS mode=${MODE}
