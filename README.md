# Cpipe  

Cpipe is a clinically focused exome sequencing pipeline developed
by the Melbourne Genomics Health Alliance. Cpipe offers an industry
standard variant calling pipeline with a suite of additional features 
needed by diagnostic laboratories added on top.

A simplified workflow is presented here. For more specific information, refer to the [documentation](docs/index.md).

## Basic Installation

To set up Cpipe, clone this repository and then run the install script:

    git clone https://github.com/MelbourneGenomics/cpipe.git
    cd cpipe
    ./install.sh
    
For more detailed instructions, have a look at the [installation documentation](docs/install.md).

## Creating your analysis batch

* Next, create the analysis directory and copy in your fastq files.
   ```bash
   mkdir -p batches/<batch_identifier>/data
   cp  <fastq_files> batches/<batch_identifier>/data
   ```
* Now, rename your fastqs to ensure they fit the following pattern:
`sampleID_<anything>_L[0-9]*_R[0-9].fastq.gz`
* Lastly, create a metadata file for your batch using:
`./cpipe batch add_batch --batch <batch_identifier> --profile ALL`

## Running the Pipeline

Now, all you need to do is run `./cpipe run --batch <batch_identifier>`
