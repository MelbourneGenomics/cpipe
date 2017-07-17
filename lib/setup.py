from setuptools import setup, find_packages
from pathlib import Path
import pkgutil

py_scripts = Path('cpipe/scripts').resolve()
sh_scripts = Path('../scripts').resolve()

entry_points = [
            f'{module} = cpipe.scripts.{module}:main'
            for _, module, _ in pkgutil.walk_packages(path=[str(py_scripts)])
            if module != 'cpipe_main'
            ] + ['cpipe = cpipe.scripts.cpipe_main:main']
scripts = [str(p) for p in sh_scripts.iterdir() if p.is_file()]

setup(
    name="Cpipe",
    version="2.5",
    packages=find_packages(),
    package_data={'cpipe.test': 'data/*'},

    # Copy the non-python scripts to the python bin directory, allowing them to be used from any directory
    scripts=scripts,

    # Create executables from a number of python modules
    entry_points={
        'console_scripts': entry_points
    },

    install_requires=[
        'requests==2.18.1',
        'python-swiftclient==3.3.0',
        'python-keystoneclient==3.12.0',
        'doit==0.30.3',
        'pymysql==0.7.11',
        'pandas==0.20.2',
        'pandas_schema==0.3.1',
        'visidata==0.91',
        'Jinja2==2.9.6',
        'argparse_prompt==0.0.2'
    ],
    license="GPL"
)
