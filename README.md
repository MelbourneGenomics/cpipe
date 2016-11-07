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
    
Ensure you follow any instructions given by the install script. For example non MGHA users must install their own
version of GATK since their license forbids us from redistributing it.
    
For more detailed instructions, have a look at the [installation documentation](docs/install.md).

## Creating your analysis batch
* Choose a batch identifier, which we will refer to as `<batch identifier>` in future commands
* Next, create the analysis directory and copy in your fastq files.
   ```bash
   mkdir -p batches/<batch identifier>/data
   cp  <fastq files> batches/<batch identifier>/data
   ```
* Now, rename your fastqs to ensure they fit the following pattern:
`sampleID_<anything>_L[0-9]*_R[0-9].fastq.gz`
* Finally, create a metadata file and config file for your batch using:
    ```bash
        ./cpipe batch add_batch --batch <batch identifier> --profile ALL --exome <target region>
    ```
    In this case, `<target region>` is the full filepath to a capture regions bed file specified by your sequencer

For more information about this stage, refer to the [batches documentation](docs/batches.md).

## Running the Pipeline

Now, all you need to do is run `./cpipe --batch <batch identifier>` run

The run command is documented in the [command documentation](docs/commands.md#run).
