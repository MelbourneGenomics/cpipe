# Dependencies
In order to run cpipe you will first need to install its dependencies, including various language runtimes and tools. These dependencies are platform-specific, so cannot be bundled with the installer or asset package.

Language dependent packages are not included since they are included in the asset bundle.

These dependencies can be automatically installed on Ubuntu 16.4 or greater using `./install-ubuntu.sh`, running the script included in this repository. If you are able to successfully do so, skip this step.

The dependencies are as follows:
* Language Runtimes 
	* Python 2.7 (including pip)
	* Perl 5.18 (including cpan)
	* R 3.2.0
	* Open JRE/JDK 1.8
* General Tools
	* git
	* make
	* gcc, g++
	* poppler-utils
	* patch
* C Libraries (probably suffixed with '-dev' or '-devel' depending on your linux distribution)
	* zlib
	* ncurses
	* openssl