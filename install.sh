#!/usr/bin/env bash

# Fail on error
set -e

# Set variables
PYTHON_VERSION="2.7.12"
ROOT=$(readlink -f $(dirname ${BASH_SOURCE}))
TEMP_PYTHON=${ROOT}/tmpdata/python
TEMP_PYBIN=${TEMP_PYTHON}/bin
PYTHON=${ROOT}/tools/python

if [[ -n $1 ]] ; then
   COMMAND=$@
else
   COMMAND=install
fi

# Load swift credentials if they exist
if [[ -f ${ROOT}/swift_credentials.sh ]] ; then
    source ${ROOT}/swift_credentials.sh
fi

# Use python-build to install python
if [[ ! -d ${PYTHON} ]]; then
	pushd tmpdata
        	git clone --depth 1 git://github.com/yyuu/pyenv.git
        	mkdir -p ${PYTHON}
        	env PREFIX=${TEMP_PYTHON} pyenv/plugins/python-build/bin/python-build ${PYTHON_VERSION} ${TEMP_PYTHON}
		rm -rf pyenv
	popd

    # Install virtualenv and create a real python installation. Activate it
    ${TEMP_PYBIN}/pip install -q virtualenv
    ${TEMP_PYBIN}/virtualenv ${PYTHON}

    # Delete the old python
    rm -rf ${ROOT}/tmpdir/*
    
fi

source ${PYTHON}/bin/activate

# Install pip dependencies
pip install -q -r requirements.txt

# Source the environment file
#source ${ROOT}/install/environment.sh

# Download assets and tools using doit
python -m doit $COMMAND
