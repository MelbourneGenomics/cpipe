# Installation

## Prerequisites
Before installing Cpipe, you need to make sure you have the following packages installed on your system:
* Java 1.8 (preferably a build later than 1.8.0_20, since Cpipe has issues with this build)
* SQLite 3 headers and library
* zlib
* bz2 lib
* lib-readline

To see a complete and up-to-date list, have a look at the `apt-get install` command in the Dockerfile

However all of these dependencies are very common and are likely on your system already. If you're using a module system, you might have to
`module load` both of these modules. If they are not available at all, they can be installed by yourself or a 
sysadmin with a command like `sudo apt-get install -y libsqlite3-dev openjdk-8-jdk`.

## Swift Credentials
If you are installing Cpipe as part of the MGHA (Melbourne Genomics Health Alliance), you have access to the Cpipe asset 
bundle, which decreases the installation time to about an hour. We intend to make this bundle accessible to everyone
in the future.
 
Download the swift credentials file in one of two ways. Either clone cloning the cpipe_util repository from github:
```bash
git clone https://github.com/MelbourneGenomics/cpipe_util
cp cpipe_util/swift_credentials.sh cpipe
```
Or download the credentials file from NECTAR (Compute → Access & Security → API Access → Download OpenStack RC File) and
place it in the cpipe root directory

## Installation
Run `./install.sh` and the full installation should take place. 

If you have the swift credentials, this should take around an hour. If you do not, this could take more than a day, depending
on the speed of your computer and its processing power

### Install Script Options
The Cpipe installer has a number of command line flags that can be used to customise the installation. You can view the 
 latest flags simply by running `./cpipe --help`. Otherwise, some are documented here:
 ```
  --help, --usage
    Print this help page to stdout
  -n, --processes <process number>
    Set the maximum number of processes to use for the install. The higher number the faster the install, but the more memory used. Defaults to the output of 'nproc --all', the number of available processing units (currently 4 on your system)
  -c, --credentials </path/to/swift_credentials.sh>
    Use the specified swift credentials file to download assets from NECTAR. Defaults to looking in the cpipe root directory
  -v, --verbose
    Print everything to stdout instead of just errors. Good for debugging an install
  -s, --no-swift
    Do a manual install instead of downloading assets from NECTAR. Strongly NOT recommended as this will potentially take days to complete
  -t, --task <taskname>
    Specify one or more tasks to run instead of a full install. Each time you use this flag it will add another task. Don't use this unless you know what you're doing
  -p, --no-pip
    Don't update pip modules. Don't use this unless you know what you're doing
 ```