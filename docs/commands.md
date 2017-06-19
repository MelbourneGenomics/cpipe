# Commands
* [Cpipe Environment](#cpipe-environment)
* [Basic Functionality](#basic-functionality)
* [Run](#run)
* [Test](#test)
* [Batch](#batch)
    * [List](#list)
    * [Create](#create)
    * [Edit](#edit)
    * [View](#view)
    * [Check](#check)
    * [Add Sample](#add-sample)
* [Visidata](#visidata)
* [Design](#design)
  
## Cpipe Environment
The Cpipe commands listed here can be run in one of two ways.

The first, simplest, way, is to use the `cpipe` script in the root of the Cpipe installation. For example,
to run a batch, you'd use:
```batch
./cpipe run MyBatch
```
In other words, you can take any command from this documentation, and replace `cpipe` with `./cpipe`.

The second way, which is useful for debugging Cpipe commands, for using additional Cpipe commands that 
aren't part of the `cpipe` command, and for checking versions of different binaries, 
is to start a new bash shell inside the Cpipe environment. You can start this shell by running the same
`./cpipe` script as above, but with no further arguments. Once you're in the shell, you can run the cpipe
commands as they're written below. For example:
```
cpipe run MyBatch
```

## Basic Functionality

The core command you will be using with is `cpipe`, which exposes the 
public interface for running and reviewing the pipeline

```
usage: cpipe [-h] {run,batch,stop,test,design,bpipe} ...

Runs core Cpipe functions

positional arguments:
  {run,batch,stop,test,design,bpipe}
    run                 Runs the analysis pipeline
    batch               Manage Cpipe batches and metadata files
    stop                Stops a currently running Cpipe batch
    test                Run the cpipe test suite
    design              Create and modify designs and their genelists
    bpipe               Pass a command to bpipe

optional arguments:
  -h, --help            show this help message and exit
```

## Run
```
usage: cpipe run [-h] [-j] batch ...

positional arguments:
  batch                The name of the batch that you want to run, relative to
                       the batches directory
  bpipe_opts           Arguments to pass to the underlying bpipe command

optional arguments:
  -h, --help           show this help message and exit
  -j, --no-java-check  Don't check that the system Java is a compatible
                       version
```
The cpipe run command is the most common use case for cpipe. Just ensure you have a sample batch in the batches
directory, run the command above, and your samples should be analysed. 

`cpipe run` can also be used to resume failed or terminated batches. Run the same command again and Cpipe should continue
where it left off.

Any arguments that you specify other than the ones above are assumed to be bpipe options. For example, if you wanted to
limit the maximum number of parallel jobs in the pipeline to 2, you could use:
```bash
cpipe run MyBatch -n 2
```
For more bpipe options, refer to the `bpipe run` documentation [here](http://bpipe-test-documentation.readthedocs.io/en/latest/Commands/run/)

## Test
```
usage: cpipe test [-h] [-j]

optional arguments:
  -h, --help           show this help message and exit
  -j, --no-java-check  Don't check that the system Java is a compatible
                       version
```

The test command runs all cpipe's unit tests and integration tests in order to ensure the pipeline
is working correctly.

## Batch
```
usage: cpipe batch [-h] [-m] {list,create,edit,view,check,add_sample} ...

positional arguments:
  {list,create,edit,view,check,add_sample}
    list                Lists the batches in the current Cpipe installation
    create              Creates a new batch, including data, metadata file and
                        configuration file
    edit                Edit the metadata file for the chosen batch
    view                View the metadatafile for the chosen batch in a human-
                        readable format
    check               Validate the metadata file for the chosen batch
    add_sample          Add a sample to an existing metadata file

optional arguments:
  -h, --help            show this help message and exit
  -m, --mgha            Use MGHA-specific validation rules
```
The batch command provides utilities for manipulating batches of samples. The 6 subcommands are listed below:

### List
```
usage: cpipe batch list [-h]

optional arguments:
  -h, --help  show this help message and exit
```
This command outputs a plaintext list of batch names. All of these batches are valid inputs for the other batch commands,
such as `edit`, `view`, and general commands like `run`.

### Create
```
usage: cpipe batch create [-h] -d DATA [DATA ...] -e EXOME [-p PROFILE]
                          [-t METADATA] [-f] [-m MODE] [-5]
                          name

positional arguments:
  name                  The name for the new batch

optional arguments:
  -h, --help            show this help message and exit
  -d DATA [DATA ...], --data DATA [DATA ...]
                        The fastq files to add to the batch
  -e EXOME, --exome EXOME
                        A bed file indicating which regions are covered by the
                        sequencing procedure
  -p PROFILE, --profile PROFILE
                        The analysis profile (gene list) to use for the
                        analysis of this batch
  -t METADATA, --metadata METADATA
                        The path to an existing metadata file to use for this
                        batch. Otherwise, generate a new metadata file.
  -f, --force           Replace an existing batch with that name, if it
                        already exists
  -m MODE, --mode MODE  Either "copy", "link" or "move": the method used to
                        put the data files into the batch directory
  -5, --no-md5          Don't check the md5 sums of fastq files
```
Creates a new batch, including adding fastqs, creating the metadata file, and creating the batch configuration file.

### Edit
```
usage: cpipe batch edit [-h] [-e EDITOR] batch

positional arguments:
  batch                 The name of the batch whose metadata file you want to
                        edit

optional arguments:
  -h, --help            show this help message and exit
  -e EDITOR, --editor EDITOR
                        The name of the executable you want to use to edit the
                        metadata file using. Defaults to visidata, included
                        with cpipe (https://github.com/saulpw/visidata)
```
`./cpipe batch show_batch <options>` 

Edits the sample metadata file as an interactive spreadsheet. You can choose a custom editor, for example `vim`, by
specifying that with the `--editor` flag.

For more information on the default visidata editor, see the [visidata](#visidata) section below.
  
### View
```
usage: cpipe batch view [-h] batch

positional arguments:
  batch       The name of the batch whose metadata file you want to view

optional arguments:
  -h, --help  show this help message and exit
```

Shows the metadata file as a spreadsheet, much like the edit command, but without any possibility of modifying the file.

### Check
```
usage: cpipe batch check [-h] batch

positional arguments:
  batch       The name of the batch whose metadata file you want to validate

optional arguments:
  -h, --help  show this help message and exit
```
Validates the given batch metadata file, using either MGHA or default validation rules (specify `--mgha` before the 
`check` command to use MGHA validation)

### Add Sample
```
usage: cpipe batch add_sample [-h] [-d DATA [DATA ...]] batch

positional arguments:
  batch                 The name of the batch to which you want to add a sample

optional arguments:
  -h, --help            show this help message and exit
  -d DATA [DATA ...], --data DATA [DATA ...]
                        The list of fastqs you want to add as samples to the
                        batch
```

Adds one or more samples to an existing batch

## Visidata
The default editor we use for editing and viewing the batch metadata is called `visidata`. You can read the official documentation on the
editor [here](https://github.com/saulpw/visidata/blob/stable/README.md). However for the sake of editing metadata, the 
only things you need to know are:
* You can move your selection around using the arrow keys
* `e` Edits the current cell. Press `Enter`/`Return` once you're done to save that cell
* `Ctrl+S` saves the current document.
* `q` quits the editor without saving

## Design
```
usage: cpipe design [-h] --profile PROFILE [--force]
                    {add_profile,show_profiles,show_genes,add_genes,remove_genes,validate,add_bed,remove_bed,show_bed}

positional arguments:
  {add_profile,show_profiles,show_genes,add_genes,remove_genes,validate,add_bed,remove_bed,show_bed}
                        command to execute

optional arguments:
  -h, --help            show this help message and exit
  --profile PROFILE     profile to update
  --force               force addition of genes
```

The `design` subcommand command provides utilities for manipulating designs. Refer to the [design](designs.md)
section for an explanation of the concept of a design.

Each command takes a design name using `--profile`, and any other input from stdin (where relevant). For example, to
add a list of genes to a design, you can run:
```
cpipe design --batch my_design < genes.txt
```
