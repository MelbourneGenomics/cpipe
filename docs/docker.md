# Docker
* [Introduction](#introduction)
* [Installation](#installation)
    * [From MGHA Docker Registry](#from-mgha-docker-registry)
    * [Building the Container Yourself](#building-the-container-yourself)
* [Running the Image](#running-the-image)
    * [Creating a Batch](#creating-a-batch)
    * [SSHing into the Container](#sshing-into-the-container)
    * [Using Docker Run](#using-docker-run)

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

## Installation
The installation step for the dockerized Cpipe involves obtaining a Cpipe *image*. Regardless of how you do this, running
containers from this image should work the same. There are two main ways of obtaining a Cpipe image:
* [Downloading from the MGHA Docker Registry](#from-mgha-docker-registry) (MGHA members only)
* [Building the container yourself](#building-the-container-yourself)

### From MGHA Docker Registry
**Please note this the docker registry is not available at the moment, please build the container yourself until we have
setup the registry**

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

That's all!

### Building the Container Yourself

In order to build the Cpipe container, follow these steps:

* 1) Clone Cpipe with:

    ```bash
    git clone https://github.com/MelbourneGenomics/cpipe --branch 2.4 --depth 1
    ```
* 2a) If you are part of MGHA, copy the swift_credentials.sh file
into the cpipe directory as explained in the [installation documentation](install.md#mgha-install).
* 2b) If you aren't part of MGHA, you'll have to manually install all the tools that we aren't able to redistribute. To
do this, follow all the instructions in the [Public Install section of the Install Documentation](install.md#public-install)
* 3) Build the container with the following command, where `<version>` is some identifier you want to tag the image with.

    ```bash
        docker build . -t cpipe:<version>
    ```
    The `<version>` tag can be the release version of Cpipe (e.g. `2.4`), or it could be the git commit hash if you think
    you will have many images from the same release (e.g. `3b592c3`)

## Running the Image

### Creating a Batch

The first step in runing the Cpipe image, as with the native install, is to create your analysis batch. However in
 Docker, this batch must contain *everything* from your filesystem that Cpipe will ever need, which includes the exome target
 file that you wouldn't normally put in the batch directory.

 To create a new batch, `cd` into a directory that
 you have write access in and that has enough space to store the results of the analysis. Once inside this directory,
 make a `batches/<batch identifier>` directory as a subdirectory:

```bash
    mkdir -p batches/<batch identifier>/data
    cp <fastq files> batches/<batch identifier>/data
    cp <target region> batches/<batch identifier>
```
For explanations of these parameters, refer to the [batch documentation](batches.md#creating-a-batch).

At this point, you have two options, you can [SSH into the Docker container](#sshing-into-the-container) and run the remaining commands as with a
native install, or you can run commands using [`docker run`](#using-docker-run)

### SSHing into the Container
Before SSHing into a docker container, make sure you are in a `screen` or `tmux` session so you are able to let the process
run in the background. This is because, if you ever type `exit` or `Ctrl+D` when inside a docker container, the container
will be killed and you'll have to start the whole process again.

To SSH into the container, run
```bash
    docker run -it --entrypoint=bash -v `pwd`/batches:/opt/cpipe/batches cpipe:<tag>
```
Here, `<tag>` is the `version` of Cpipe you pulled from the registry, or the `tag` you gave to cpipe when running `docker build`.

This should place you in the cpipe directory of a fully configured Cpipe installation.

However, before you can start the analysis, you'll first have to finish creating your metadata and configuration files for your batch. To do so, run:
```bash
  ./cpipe batch add_batch --batch <batch identifier> --profile <profile name> --exome <target region>
```

All done! You can now run any of the commands listed in the [commands documentation](commands.md), including `./cpipe run`,
the main analysis command.

### Using Docker Run

To create the metadata and configuration file from outside the container, run the following command.
Here, the `<tag>` is the `version` of Cpipe you pulled from the registry, or the `tag` you gave to cpipe when running `docker build`.
```bash
    docker -v `pwd`/batches:/opt/cpipe/batches run cpipe:<tag> batch add_batch --batch <batch identifier> --profile ALL --exome <target region>
```

What the `-v` flag
does is mount your `batches` directory into the container at `/opt/cpipe/batches` as though it were a drive. The trick here is that when you're using
filepaths inside the container, you'll either have to use relative paths, which will be relative to the cpipe installation
directory (e.g. `./batches/MY_BATCH`) will be the path to the `MY_BATCH` batch. Alternatively you can use absolute paths,
for example, the absolute path to `MY_BATCH` would be `/opt/cpipe/batches/MY_BATCH`.

Because of this, in the command above, you'll have to specify the location of the <target region>` file from inside the container. For this we
recommend using the relative path to the file, for example `--exome batches/MY_BATCH/MY_REGIONS.bed`

At this point, you can start using the docker image in the exact same way as the `./cpipe` command. The only difference is
that you'll need to replace the `./cpipe` command with the same `docker run` command we used above. For example, to run
the analysis, you'd type:
```bash
    docker run  -v `pwd`/batches:/opt/cpipe/batches cpipe:<tag> run
```

If this incantation is a bit too much, you might consider aliasing a cpipe command to this large command:
```bash
alias cpipe=docker run  -v `pwd`/batches:/opt/cpipe/batches cpipe:<tag>
```

Also recall how filepaths work inside the container as explained above.

You can now refer to the [command documentation](commands.md), remembering that you should use the above command instead
of `./cpipe`.
