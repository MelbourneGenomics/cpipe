# Installation

* [Prerequisites](#prerequisites)
  * [Runtimes](#runtimes)
  * [Libraries](#libraries)
  * [Build Tools](#build-tools)
  * [Utilities](#utilities)
* [Swift Credentials](#swift-credentials)
* [Installation](#installation-1)
  * [Install Script Options](#install-script-options)
* [Set Installation Name (optional)](#set-installation-name-optional)

## Prerequisites
Before installing Cpipe, you need to make sure you have the following packages installed on your system. Most of these
are very common packages and probably came pre-installed on your machine. However this is an exhaustive list:

### Runtimes
* Java 1.8 (preferably a build later than 1.8.0_20, since Cpipe has issues with this build)

### Libraries
* SQLite 3
* zlib
* bzip2
* readline
* ssl
* xorg
* curl
* ncurses

### Build Tools
* gcc
* gfortran
* make

### Utilities
* curl
* git 

To see a complete and up-to-date list, have a look at the `apt-get install` command in the Dockerfile

If you're using a module system, you might have to
`module load` some of these modules. If they are not available at all, they can be installed by yourself or a 
sysadmin with a command like `sudo apt-get install -y libsqlite3-dev` etc..

## Installation Type
There are two official ways to install Cpipe: as a member of the MGHA (Melbourne Genomics Health Alliance), or as a 
member of the public. Please read the section relevant to you

### MGHA Install
If you are installing Cpipe as part of the MGHA, you are covered under our license for various software, so don't need
to obtain your own copies of tools like GATK.

In order to make use of this, you'll need to obtain the MGHA Swift credentials, and place them in the cpipe directory
as a file named `swift_credentials.sh`. These credentials can be obtained in one of two ways:
* Gaining read access to the cpipe_util GitHub repository (contact the dev team on help@melbournegenomics.org.au for 
access) and running the following
 commands from the parent directory **above** where Cpipe is installed (the directory that contains a cpipe subdirectory):
   ```bash
       git clone https://github.com/MelbourneGenomics/cpipe_util
       cp cpipe_util/swift_credentials.sh cpipe
   ```
   If you cloned cpipe using a different name, copy into the directory you created instead of `cpipe`
* Having access to the MGHA NeCTAR project and downloading the credentials file from:
 Compute (in the sidebar) → Access & Security → API Access → Download OpenStack RC File.
 You'll then have to rename the file `Melbourne_Genomics_Health_Alliance-openrc.sh` to `swift_credentials` and place it
 in the Cpipe installation directory.
 
### Public Install
If you *aren't* a member of the MGHA, you don't need to perform any steps prior to running the install script. However
the script itself will prompt you to download some extra tools that we can't legally distribute to you. Make sure you
read all of the install script's output to see this.

Currently, the only tool you'll need to install yourself is GATK. You can obtain the jar file from the 
[GATK Download Page](https://software.broadinstitute.org/gatk/download/). You'll then need to place it in 
`cpipe/tools/java_libs`.

## Installation
Run `./install.sh` and the full installation should take place. 

Depending on the speed of your connection, and the number of simultaneous processes your machine is able to run,
 this process could take as little as 30 minutes, but may well take longer as it involves download several gigabytes of data. 

### Install Script Options
The Cpipe installer has a number of command line flags that can be used to customise the installation. You can view the 
 latest flags simply by running `./cpipe --help`. Otherwise, some are documented here:
 ```
  --help, --usage
    Print this help page to stdout
  -n, --processes <process number>
    Set the maximum number of processes to use for the install. The higher number the faster the install, but the more memory used. Defaults to the output of 'nproc --all', the number of available processing units (currently 4 on your system)
  -c, --credentials </path/to/swift_credentials.sh>
    Use the specified swift credentials file to download assets from NeCTAR. Defaults to looking in the cpipe root directory
  -v, --verbose
    Print everything to stdout instead of just errors. Good for debugging an install
  -s, --no-swift
    Do a manual install instead of downloading assets from NeCTAR. Strongly NOT recommended as this will potentially take days to complete
  -t, --task <taskname>
    Specify one or more tasks to run instead of a full install. Each time you use this flag it will add another task. Don't use this unless you know what you're doing
  -p, --no-pip
    Don't update pip modules. Don't use this unless you know what you're doing
 ```
 
## Set Installation Name (optional)
Once the installation is completed, you can set the name of your installation, which will become a prefix on all of the pipeline output files,
 ensuring they are unique to your installation. 
 
 To do this, create a file named `pipeline_id` in the root of your Cpipe installation, and inside the file write something like
 `<sitename>_<version>` where `sitename` is the name of your system/site, and version is the release version of Cpipe you
 are using. For example, `vlsci_2.3`. You can do this in with the following command:
 ```bash
 echo "vlsci_2.3" > pipeline_id
 ```
