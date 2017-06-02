# Docker
* [Introduction](#introduction)
* [Obtaining the Container](#obtaining-the-container)
	* [From MGHA Docker Registry](#from-mgha-docker-registry)
	* [Building the Container Yourself](#building-the-container-yourself)
* [Running the Image](#running-the-image)

## Introduction
Docker is a containerization system that will allow you to install and run Cpipe without any dependencies or installation.

For this reason, Docker is the recommended method for installing Cpipe on small computing systems (e.g. cloud computing nodes).
However since Cpipe cannot easily interface with a HPC queueing system (slurm, torque etc.) from within Docker, we
currently recommend that you perform a native install on these systems.

Docker consists of *images* which are self contained filesystems that can be easily distributed and shared. Once you have
an image, you can run it, which will create a *container*, which is a mini operating system you can issue commands to
and even SSH into. For more information on how Docker works, refer to the [Docker website](https://www.docker.com/what-docker)

If you are planning on running Cpipe in a docker container, you can follow these instructions instead of those in the
[README](../README.md).

## Obtaining the Container
The installation step for the dockerized Cpipe involves obtaining a Cpipe *image*. Regardless of how you do this, running
containers from this image should work the same. There are two main ways of obtaining a Cpipe image:
* [Downloading from the MGHA Docker Registry](#from-mgha-docker-registry) (For MGHA members)
* [Building the container yourself](#building-the-container-yourself) (For everyone else)

### From MGHA Docker Registry
The easiest way to obtain a Cpipe image is by logging onto the Cpipe docker registry and downloading to the image. However
since the images contain licensed software, we can unfortunately only provide these images to MGHA members. If you are
a member of the MGHA and would like to obtain docker registry credentials, please send an email to help@melbournegenomics.org.au.

Once you have the credentials, you'll first need to login to our registry. Insert the credentials as prompted.
```bash
docker login https://docker.melbournegenomics.org
```

Now all you need to do is run the following command, where `<version>` is the version of Cpipe you would like to install
```bash
docker pull docker.melbournegenomics.org/cpipe:<version>
```
The current versions we have available on the docker registry are the following. Any of these version numbers can be put 
after the colon in the `docker pull` command above:
* 2.4.0
* 2.4.1

### Building the Container Yourself

In order to build the Cpipe container, follow these steps:

1. Clone Cpipe with:

    ```bash
    git clone https://github.com/MelbourneGenomics/cpipe --branch 2.4 --depth 1
    ```
2.
    1. If you are part of MGHA, copy the swift_credentials.sh file into the cpipe directory as explained in the [installation documentation](install.md#mgha-install).
    2. If you aren't part of MGHA, you'll have to manually install all the tools that we aren't able to redistribute. To do this, follow all the instructions in the [Public Install section of the Install Documentation](install.md#public-install)
3. `cd` into the cpipe directory and build the container with the following commands,
 where `<version>` is some identifier you want to tag the image with.
    ```bash
    cd cpipe
    docker build . -t cpipe:<version>
    ```
    The `<version>` tag can be the release version of Cpipe (e.g. `2.4`), or it could be the git commit hash if you think
    you will have many images from the same release (e.g. `3b592c3`)

## Running the Image

1. **Gathering the Input Files**: The first step in running the Cpipe image is to gather all the files you'll be using as inputs into one directory. In 
	that directory you should put all of your fastq files, and the target region BED file (refer to the `EXOME_TARGET` variable
	in [the documentation](configuration.md#options)). It may be tempting to symlink (`ln -s`) these files into the directory
	but symlinks won't work in a Docker container, so copy them instead. Your input directory should look something like
	this:
	```
	exons.bed
	NA12878_CARDIACM_MUTATED_L001_R1.fastq.gz
	NA12878_CARDIACM_MUTATED_L001_R2.fastq.gz
	```
2.  **Choosing an Output Directory**: Next, you'll need to choose a directory that you have write access in and that has enough space to store the results of 
	the analysis. We'll refer to this as the output directory
	
3.  **Starting the Container**: Now that you have all your files, you can start the container. Before you do this, however,
 make sure you're in a `screen` or `tmux` session so that the container won't be killed when you disconnect from the server
  hosting Docker. Run the following command to start the container:
	```bash
	  docker run -it -v /path/to/input/:/input -v /path/to/output:/opt/cpipe/batches/001 cpipe:<tag>
	```
	This should put you in an interactive bash shell, inside the container, similar to if you'd run the 
	[environment shell](batches.md#cpipe_environment).

4. **Creating a Batch**:Now that you're in the container, you can create a batch with the ordinary `cpipe batch` command. For example:
	```bash
	cpipe batch create --data /input/*.fastq.gz --exome /input/*.bed
	```
	
5. **Starting the Analysis**: Now that the batch has been created, all you need to do is start the analysis:
	```bash
	cpipe run 001
	```
6. **Getting the Results**: Once the analysis has finished successfully, you will find the results inside the output
	directory that you mounted into the container. If you are happy with these results, you can kill the container by
	typing `exit` in the docker shell.

