#!/usr/bin/env bash

# Fail on error
set -e

# Set variables
PYTHON_VERSION="2.7.12"
ROOT=$(readlink -f $(dirname ${BASH_SOURCE}))
PYTHON=${ROOT}/tools/python
PYBIN=${PYTHON}/bin

# Load swift credentials if they exist
if [[ -f ${ROOT}/swift_credentials.sh ]] ; then
    source ${ROOT}/swift_credentials.sh
fi

# Use python-build to install python
if [[ ! -d ${PYTHON} ]]; then
	pushd tmpdata
        	git clone --depth 1 git://github.com/yyuu/pyenv.git
        	mkdir -p ${PYTHON}
        	pyenv/plugins/python-build/bin/python-build ${PYTHON_VERSION} ${PYTHON}
		rm -rf pyenv
	popd
fi

# Install python dependencies
${PYBIN}/pip install -r requirements.txt

# Download assets and tools using doit
${PYBIN}/python -m doit install

# Export the python location so this can be source'd
export PATH=${PYBIN}:$PATH