#!/usr/bin/env bash

# Fail on error
set -e

# Set variables
PYTHON_VERSION="2.7.12"
ROOT=$(readlink -f $(dirname ${BASH_SOURCE}))
PYTHON=${ROOT}/tools/python

if [[ ! -f ${PYTHON}/python ]]; then
	
	# Install Python
	mkdir -p ${PYTHON}
	curl https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz | tar -xzf - --strip-components=1 -C ${PYTHON}
	pushd ${PYTHON}
    		./configure
    		make
	popd
	
fi

# Install pip
${PYTHON}/python -m ensurepip

# Install python dependencies
${PYTHON}/pip install -r requirements.txt

# Download assets and tools using doit
${PYTHON}/python -m doit install

# Export the python location so this can be source'd
export PATH=${PYTHON}:$PATH

