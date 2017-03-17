import os
from setuptools import setup, find_packages
from pathlib import Path
import pkgutil
import cpipe.scripts

# root = pathlib.Path(__file__).parent.parent.resolve()
scripts = Path('scripts')

entry_points = []
for importer, module, ispkg in pkgutil.walk_packages(path=cpipe.scripts.__path__):
    entry_points.append(f'{module} = cpipe.scripts.{module}:main')

setup(
    name="Cpipe",
    version="2.5",
    packages=find_packages(),
    package_data={'cpipe.test': 'data/*'},

    # Copy the non-python scripts to the python bin directory, allowing them to be used from any directory
    scripts=[str(p) for p in scripts.iterdir() if p.is_file()],

    # Create executables from a number of python modules
    entry_points={
        'console_scripts': entry_points
    },

    install_requires=[
        'requests==2.11.0',
        'python-swiftclient',
        'python-keystoneclient',
        'doit>=0.30.0',
        'pymysql',
        'pandas',
        'pandas_schema',
        'visidata>=0.37'
    ],
    license="GPL"
)
