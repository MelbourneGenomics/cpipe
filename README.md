# Cpipe  

Cpipe is a clinically focused exome sequencing pipeline developed
by the Melbourne Genomics Health Alliance. Cpipe offers an industry
standard variant calling pipeline with a suite of additional features 
needed by diagnostic laboratories added on top.

A simplified workflow is presented here. For more specific information, refer to the [documentation](docs/index.md).

## Installation

To set up Cpipe, clone this repository and then run the install script:

    git clone https://github.com/MelbourneGenomics/cpipe.git
    cd cpipe
    ./install.sh
    
Ensure you follow any instructions given by the install script. For example, non MGHA users must install their own
version of GATK since their license forbids us from redistributing it.
    
For more detailed instructions, have a look at the [installation documentation](docs/install.md).

## Activate the Cpipe environment
```bash
./environment.sh
```
This will start a new shell with various Cpipe-specific variables set. Just press `Ctrl+D` or type `exit` to close this
shell at any time.

## Creating your analysis batch
* First, rename your fastqs to ensure they fit the following pattern:
`sampleID_<anything>_L[0-9]*_R[0-9].fastq.gz`
* Tell Cpipe to create a new batch using your fastq files
   ```bash
   manage_batch create MyBatch --data path/to/samples/*.fastq.gz --exome path/to/exons.bed
   ```
  * `MyBatch` is the identifier 
  * `path/to/exons.bed` is the full filepath to a capture regions bed file specified by your sequencer
  * `path/to/samples/*.fastq.gz` is the full filepath to the sample fastqs you want to put in your batch

For more information about this stage, refer to the [batches documentation](docs/batches.md).

## Running the Pipeline

Now, all you need to do is run `cpipe --batch MyBatch run`

The run command is documented in the [command documentation](docs/commands.md#run).
