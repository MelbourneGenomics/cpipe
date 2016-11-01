# Running Cpipe
* [Basic Functionality](#basic-functionality)
* [Run](#run)
* [Test](#test)
* [Batch](#batch)
  * [Add Batch](#add-batch)
  * [Show Batch](#show-batch)
  * [Show Batches](#show-batches)
  * [Add Sample](#add-sample)
* [Genelist](#genelist)
* [Metadata](#metadata)
  * [Update](#update)
  * [Validate](#validate)

## Basic Functionality

The basic Cpipe entrypoint is the cpipe script, which exposes all the public interface for running and reviewing the 
pipeline

Usage: `./cpipe <CPIPE OPTIONS> COMMAND <COMMAND OPTIONS>`

Commands (type --help after any command for more details):
* [run](#run): Runs the analysis pipeline
* [test](#test): Runs the pipeline tests
* [batch](#batch): Creates and modifies analysis batches
* [genelist](#genelist): Creates and modifies genelists
* [metadata](#metadata): Creates and modifies sample metadata files

Keyword Arguments:
* `-b, --batch <batch name>`: Specify a batch (a subdirectory inside batches) to use for the run and bpipe commands. Defaults to a batch named 'batch'
* `--help, --usage`: Prints this help page
* `-j, --no-java-check`: Disables the java version check. Only do this if you know what you're doing

## Run
The cpipe run command is the most common use case for cpipe. Just ensure you have a sample batch in the batches directory 
and run `./cpipe run --batch <batch-name>` and your samples should be analysed. 

You can also specify bpipe options (like limiting the system resources) with the `--bpipe-options <OPTS>` flag

## Test
The test command is also very useful. It runs all cpipe's unit tests and integration tests in order to ensure the pipeline
is working correctly. Test takes no additional options.

## Batch
The batch command provides utilities for manipulating batches of samples. The batch subcommands are:
* add_batch
* show_batches
* show_batch
* add_sample

### Add Batch
`./cpipe batch add_batch <options>`

Creates a new sample metadata file for a new batch

Options:
  * `--batch BATCH` (required) Specifies the batch directory in which to create a metadata file
  * `--profile PROFILE` (optional) Specifies the analysis profile to use. Defaults to ALL, the recommended profile
  * `--exome EXOME` (optional) Specifies a bed file to use as the capture region
  * `--data [DATA [DATA ...]]` (optional) specifies the fastq files in the batch. Defaults to all files in the data 
  subdirectory, which is the recommended structure for a batch.
  * `--force` (optional) Overrides an existing samples.txt

### Show Batch
`./cpipe batch show_batch <options>` 

Prints information about a single batch

Options:
  * `--batch BATCH` (required) Specifies the batch to print information about
  
### Show Batches
`./cpipe batch show_batches <options>`

Lists all the existing batches in this cpipe installation

Options:
  * `--batch BATCH` (required) Specifies the batch to print information about

### Add Sample
`./cpipe batch add_sample <options>`

Adds one or more samples to an existing batch

  * `--batch BATCH` (required) Specifies the batch to add a sample to
  * `--profile PROFILE` (optional) Specifies the analysis profile to use. Defaults to ALL, the recommended profile
  * `--data [DATA [DATA ...]]` (optional) specifies the fastq files to add to the in the batch. Defaults to all files in the data 
  subdirectory, which is the recommended structure for a batch.

## Genelist

The genelist command provides utilities for manipulating designs (see the terminology section for an explanation 
of this term).

The genelist command has the following subcommands:

* `add_profile` 
* `show_profiles` 
* `show_genes` 
* `add_genes < genes.txt [--force]` 
* `remove_genes < genes.txt [--force]` 
* `validate` 
* `add_bed < genes.bed [--force]` 
* `remove_bed < genes.bed` 
* `show_bed` 

All commands have the following arguments:
 * `--profile PROFILE` (mandatory): Select a profile/design to update 
 * `-h, --help`: show the help message and exit

## Metadata

The metadata command provides utilities for manipulating metadata files, which define the options for the analysis of a 
given batch. The metadata command has two subcommands: `check` for checking the validity of an existing metadata file,
and `update`, for modifying an existing metadata file. Refer to the add_batch subcommand of the [`batch`](#batch) command 
for a way of creating a new metadata file.

### Update
Usage: `./cpipe metadata check < samples.txt --fields [field] [field...]`
Prints out the values of the metadata file. Only prints the fields specified if `--fields` is present

### Validate
Usage: `./cpipe metadata validate`
Updates a given field of the sample metadata file.

Arguments:
* `--sample_id <sample>`: Specifies the sample ID of the sample to update
* `--name <name>`: specifies the name of the field to update
* `--value <value>`: specifies the new value to assign the field
* `--target <file>`: specifies the sample metadata file to update
