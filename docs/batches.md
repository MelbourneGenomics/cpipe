# Batches

* [Introduction](#introduction)
* [Manipulating Batches](#manipulating-batches)
  * [Creating a Batch](#creating-a-batch)
* [Files](#files)
  * [Data Directory](#data-directory)
  * [Sample Metadata](#sample-metadata)
* [config.batch.groovy](#configbatchgroovy)

## Introduction
In Cpipe, a batch is a group of samples to be analysed at the same time. In the filesystem, a batch is a directory inside
 `cpipe/batches`. For example, a batch named `batch_001` would mean creating a `cpipe/batches/batch_001`gi directory.
 
## Manipulating Batches
### Creating a Batch
Once you have your fastq files, follow these steps to create a new analysis batch:
* Create the batch directory and copy in the fastq data:
    ```bash
    mkdir -p batches/<batch identifier>/data
    cp <fastq files> batches/<batch identifier>/data
    ```
* Create the metadata file using:
    ```bash
      ./cpipe batch add_batch --batch <batch identifier> --profile <profile name> --exome <target region>
    ```
    * `<target region>` is the full filepath to a capture regions bed file specified by your sequencer. For example, for
    Illumina sequencing prepared with a Nextera DNA Library Preparation Kit, this file can be downloaded from 
    [the Illumina website](http://support.illumina.com/downloads/nextera-rapid-capture-exome-v1-2-product-files.html).
    * `<profile name>` is the name of the design (subdirectory within `cpipe/designs`) to use for the analysis. By default
    this is the ALL profile, which we recommend using. See the [designs documentation](designs.md) for more information.

  For more information on the `add_batch` command, refer to its [documentation](commands.md#add-batch)
  
### Creating a Batch for Re-analysis
If you intend to re-analyse an existing batch from an old version of Cpipe, first create the new batch directory, then
copy the data directory, sample metadata file, and config file into this new batch. You can use this script to assist 
you.

First set the `OLD_BATCH` and `NEW_BATCH` variables to the path to the relevant batch directory (e.g. `OLD_BATCH=batches/001`, 
`NEW_BATCH=/path/to/old/cpipe/batches/001`
and then run this script:
```bash
mkdir -p $NEW_BATCH
for file in data/ samples.txt target_regions.txt config.batch.groovy ; do
    if [[ -e $OLD_BATCH/$file ]] ; then
        cp -R $OLD_BATCH/$file $NEW_BATCH
    fi
done
```

Note that you don't need to run the `add_batch` script as you do when creating a new batch, because this creates the
 metadata file and configuration file which should already exist if you are performing reanalysis.
  
### Adding More Samples

To add more samples to an existing batch, use the `./cpipe batch add_samples` command. Refer to 
[its documentation](./commands.md#add-sample)

### Viewing Batch Information

Cpipe provides two utility commands for viewing batch information:
* `./cpipe batch show_batches` lists all the batches in the current installation. Refer to 
[its documentation](commands.md#show-batches) for more information
* `./cpipe batch show_batch --batch <batch name>` will list information about an existing batch. Refer to
[its documentation](commands.md#show-batch) for more information.

## Files
Inside the batch directory (`batches/<batch name>`) are three fundamental elements. All of these must be present and 
are not optional. Here they are covered in more detail.
 * The `data` subdirectory
 * The sample metadata file
 * A `config.batch.groovy` configuration file
  
### Data Directory 
The data directory is a directory named `data` that will hold all of the fastq samples for this batch. 
Each sample in this directory must fit the pattern `<sample id>_<anything>_<lane number>_<read number>.fastq.gz`, for 
example, `00NA12877_Coriell_000_TGx140395_TL140776_L001_R1.fastq.gz`. The components to this filename are as follows
* `sample id` is any unique name for the sample, e.g. `NA12878` in our example
* `anything` of course can be anything produced by the sequencer, in this case `Coriell_000_TGx140395_TL140776` 
* `lane number` must be the letter L followed by the lane number, in this case L001
* `read number` must be the letter R followed by either 1 or 2, the read number, which is `R1` in this example

### Sample Metadata
The sample metadata is the `samples.txt` located in the batch directory (batches/<batch name>) of the batch you're 
running. The file is a TSV (tab-separated text file), where each line corresponds to an input sample. Here you can set
various options about each sample

|Position | Name | Notes | Allowed Values | Availability|
|---|---|---|---|---|
|1 | Batch | Identifier for this batch (for Melbourne Genomics this is a 3-digit ID) | Unrestricted | 2.2+|
|2 | Sample ID | Identifier for this sample (for Melbourne Genomics this is a 9-digit ID) | Unrestricted | 2.2+|
|3 | DNA Tube ID |  	 | Unrestricted | 2.2+|
|4 | Sex |  	 | Male, Female, Unknown, other | 2.2+|
|5 | DNA Concentration | ng/uL | Numeric [0-9.] | 2.2+|
|6 | DNA Volume | uL | Unrestricted | 2.2+|
|7 | DNA Quantity | ng | Numeric [0-9.] | 2.2+|
|8 | DNA Quality |  	 | Numeric [0-9.] | 2.2+|
|9 | DNA Date | Date of extraction | Comma separated list of dates of the format  yyyymmdd | 2.2+|
|10 | Cohort | Name of target region to be analysed for the patient | Unrestricted | 2.2+|
|11 | Sample Type |  	 | Normal, Tumour | 2.2+|
|12 | Fastq files |  	 | Comma separated list of FASTQ files found in the data directory. | 2.2+|
|13 | Prioritised genes | The categories are used to prioritise variants from these genes in the gene priority column of the pipeline output. | Comma separated gene list, space separating priorities, e.g. 3:GABRD,KCNAB2,ALG6 4:CASQ2,HAX1,CHRNB2,KCNJ10 | 2.2+ |
|14 | Consanguinity |  	 | No, Yes, Suspected, Unknown| 2.2+ |
|15 | Variants file | Known variants for the disease (not implemented) | Unrestricted |2.2+|
|16 | Pedigree file | PED specification for trios (see Trio Analysis below) | Unrestricted|2.2+|
|17 | Ethnicity | For filtering on specific variants (not implemented) | Unknown, European, African, Asian|2.2+ |
|18 | Variant call group | Comma separated list of samples to call as a group (not implemented) | Unrestricted| 2.2+ |
|19 | Capture date | Exome capture date | Comma separated list of dates of the format  yyyymmdd|2.2+
|20 | Sequencing Date | Date of sequencing. | Comma separated list of dates of the format  yyyymmdd|2.2+ |
|21 | Mean Coverage | Total on-target aligned mean coverage, post duplicate removal as obtained by the sequencing laboratory. | Numeric [0-9.]|2.2+ |
|22 | Duplicate % | Percentage of detected duplicates removed before calculating mean on-target coverage | 2.2+ |
|23 | Machine ID | Provided by the sequencing laboratory | Comma separated list | 2.2+|
|24 | DNA Extraction Lab | Lab that extracted the DNA | Unrestricted | 2.2+ |
|25 | Sequencing Lab | Lab that sequenced the sample | Unrestricted | 2.2+ |
|26 | Exome capture | Capture kit used (e.g. Nextera Rapid Capture Exome 1.2) | Unrestricted | 2.2+ |
|27 | Library preparation | Additional library preparation | Unrestricted | 2.2+ |
|28 | Barcode pool size | Size of barcode adapters | Unrestricted | 2.2+|
|29 | Read type | The type of read (e.g. 100PE)  | Unrestricted | 2.2+ |
|30 | Machine type | Sequencing machine (e.g. HiSeq 2500) | Unrestricted | 2.2+ |
|31 | Sequencing chemistry | e.g. TruSeq Rapid XPS | Unrestricted | 2.2+ |
|32 | Sequencing software | e.g. Illumina RTA 1.18.64 | Unrestricted | 2.2+ |
|33 | Demultiplex software | e.g. bcl2fastq 1.8.4 | Unrestricted | 2.2+|
|34 | Hospital centre | Origin of the patient sample | Unrestricted | 2.2+ |
|35 | Sequencing contact | Where sequencing alerts should be sent | Unrestricted | 2.2+ | 
|36 | Pipeline contact | Where pipeline result alerts should be sent | Unrestricted | 2.2+ | 
|37 | Notes | Additional notes or relevant information about the sequencing | Unrestricted | 2.2+ |
|38 | Pipeline notes | Additional notes relevant to the operation of the pipeline | Unrestricted | 2.3 |
|39 | Analysis type | Currently unused | Unrestricted | 2.3 |

### config.batch.groovy
The final file that must be included in the batch directory is the configuration file. Refer to the [configuration](configuration.md)
section for more details.

Note that, while most of the fields are optional in this file, the config.batch.groovy **must** contain an `EXOME_TARGET`
field. This should be set automatically by running `./cpipe batch add_batch` as documented in the [Creating a Batch](#creating-a-batch)
section.
