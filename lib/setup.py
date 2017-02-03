import os
from setuptools import setup, find_packages
from pathlib import Path

# root = pathlib.Path(__file__).parent.parent.resolve()
scripts = Path('cpipe') / 'scripts'

setup(
    name="Cpipe",
    version="2.5",
    packages=find_packages(),

    # Copy the non-python scripts to the python bin directory, allowing them to be used from any directory
    # scripts=[str(p) for p in scripts.iterdir() if p.is_file()],

    # Create executables from a number of python modules
    entry_points={
        'console_scripts': [
            'cpipe = cpipe.scripts.cpipe',
        ]
    },
    install_requires=[
        'python-swiftclient',
        'python-keystoneclient',
        'doit==0.30.0',
        'pymysql',
        'pandas',
        'pandas_schema',
        'visidata==0.37'
    ],
    license="GPL"
)
