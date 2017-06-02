# Development
## Project Structure
The main directories within the source directory are:

* **`batches`**: Contains the input and output files from analysis. This is mostly the domain of users and not developers. Refer to the [batches documentation](batches.md) for more information.
* **`designs`**: Contains the design directories, which are also for user modification, except for the built in ALL and CARDIOM designs, which come with Cpipe. In addition, the genelists directory contains shared region information that might need occasional updating. See the [designs documentation](designs.md) for more information.
* **`docs`**: Contains documentation in markdown format.
* **`lib`**: Contains the Python source code for Cpipe, as a python package. Python-specific tests and scripts are also kept here
* **`pipeline`**: Contains the Groovy source code, namely the core of the pipeline, including tests, and the main pipeline entry point, pipeline.groovy
* **`scripts`**: Interpreted scripts written in languages other than python, that are used by the pipeline or as utility scripts. Includes bash scripts, groovy scripts etc, with the relevant shebang line. Note that these scripts should **not** have a file extension, because this means the binary created from these will also have that extension, which does not meet the unix standard practise. These scripts are copied into the local bin directory by the install process.
* **`tasks`**: Contains tasks written in Python using the [doit](http://pydoit.org/contents.html) automation framework. Currently only used for installation.
* **`tmpdata`**: Used instead of /tmp for the install process, as there is likely more storage in the drive Cpipe is being installed into, rather than the root filesystem.
* **`tools`**: Will contain the tools installed during the installation process.
* **`data`**: Will contain the genomic reference data installed during the installation process.
