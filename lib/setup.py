from setuptools import setup, find_packages
from pathlib import Path
import pkgutil

py_scripts = Path('cpipe/scripts').resolve()
sh_scripts = Path('scripts').resolve()

entry_points = [
            f'{module} = cpipe.scripts.{module}:main'
            for _, module, _ in pkgutil.walk_packages(path=[str(py_scripts)])
            if module != 'cpipe_main'
            ] + ['cpipe = cpipe.scripts.cpipe_main:main']
scripts = [str(p) for p in sh_scripts.iterdir() if p.is_file()]

print('\n'.join(scripts))
print('\n'.join(entry_points))

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
