#!/usr/bin/env bash

# Fail on error
set -e

# Set variables
PYTHON_VERSION="2.7.12"
ROOT=$(readlink -f $(dirname ${BASH_SOURCE}))
PYTHON=${ROOT}/tools/python

echo $BASH_SOURCE $ROOT

# Install Python
mkdir ${PYTHON}
curl https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz | tar -xzf --strip-components=1 -C ${PYTHON}
pushd ${PYTHON}
    make
popd
export PATH=${PYTHON}/bin:$PATH

# Install python dependencies
pip install -r requirements.txt

# Download assets and tools using doit
doit install
