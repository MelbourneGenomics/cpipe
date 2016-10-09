#!/usr/bin/env bash

# Fail on error
set -e

# Set variables
PYTHON_VERSION="2.7.12"
ROOT=$(readlink -f $(dirname ${BASH_SOURCE}))
TEMP_PYTHON=${ROOT}/tmpdata/python
TEMP_PYBIN=${TEMP_PYTHON}/bin
PYTHON=${ROOT}/tools/python

# Load swift credentials if they exist
if [[ -f ${ROOT}/swift_credentials.sh ]] ; then
    source ${ROOT}/swift_credentials.sh
fi

# Use python-build to install python
if [[ ! -d ${PYTHON} ]]; then
	pushd tmpdata
        	git clone --depth 1 git://github.com/yyuu/pyenv.git
        	mkdir -p ${PYTHON}
        	pyenv/plugins/python-build/bin/python-build ${PYTHON_VERSION} ${TEMP_PYTHON}
		rm -rf pyenv
	popd
fi

# Install virtualenv and create a real python installation. Activate it
${TEMP_PYBIN}/pip install virtualenv
${TEMP_PYBIN}/virtualenv ${PYTHON}
source ${PYTHON}/bin/activate

# Source the environment file
#source ${ROOT}/install/environment.sh

# Delete the old python
rm -rf ${ROOT}/tmpdir/*

# Install pip dependencies
pip install -r requirements.txt

# Download assets and tools using doit
python -m doit install
