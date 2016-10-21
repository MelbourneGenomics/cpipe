# Configuration

The global configuration of Cpipe can be modified in 2 different ways:

* [By modifying the `config.groovy` or `config.batch.groovy` files](#config.groovy)
* [By modifying `bpipe.config`](#bpipe.config)

We will explore these options one at a time.

## config.groovy

The config.groovy files are the main configuration files for cpipe. 

The main `config.groovy` file, located in the pipeline
directory, has settings that affect every batch run with this installation. For this reason, we strongly recommend that 
you avoid changing this file, because if you change this file between two sets of analysis, it will be unclear with which
settings the first was run. The only possible time you might want to edit this file is to change one of the tools 
directories, but since Cpipe bundles all of its own tools, this isn't recommended either.

Instead, we suggest you create a `config.batch.groovy` file in your batch directory (batches/<batch name>) and alter
configuration there. 

### Syntax

* Each option must be provided as though declaring a shell variable, namely in capital letters with
no space around the equals sign
* Comments are also allowed in these files, but they must be using C-style slashes rather than bash hash symbols
* For example, these lines are valid for the config.groovy file:
    ```
    // Whether to fail analysis if FASTQC produces warnings
    CHECK_FASTQC_FAILURES=false
    ```
    
### Options

The most up-to-date list of config options can be found simply by looking through the pipeline/config.groovy file in the
source code. However we will also document those variables here.

| Option | Description | Editable |
| --- | --- | --- |
| BASE | The base of everything - set this to the absolute path of the root of the pipeline distribution (most likely, parent folder of the folder this file is in) | No |
| DATA | The base directory for reference data. This parameter isn't used directly by the pipeline, but most reference file locations are defined in terms of this variable | No |
| TMPDIR | The directory in which to store temporary files | Yes |
| EMAILS | An email address which will be notified about failures | Yes |
| REF | The location of the genome reference FASTA file | No |
| DBSNP | The location of the DBSNP VCF file | No |
| GOLD_STANDARD_INDELS | The location of the VCF file containing known indels | No |
| EXOME_TARGET | The location of a BED file specifying the capture region of your sequencer. For example, the Nextera capture region files can be downloaded from the Illumina website: http://support.illumina.com/sequencing/sequencing_kits/nextera-rapid-capture-exome-kit/downloads.html | Yes |
| SAMPLE_ID_MASK | The variant database requires every sample id to map to a different individual. If multiple sample ids can map to the same individual (for example, you repeat sequencing on a sample, etc.), then you can mask out the part of the sample id that is not unique with a regular expression here. For example, Melbourne Genomics uses a 9 digit sample (study) id, but only the first 7 digits are unique for an individual. We can use a mask of ".{7}" to indicate that only the first 8 digits of the study id should be used in the variant database. | Yes |
| OUT_OF_COHORT_VARIANT_COUNT_FILTER | Filter out variants observed more than this many times (ie: if this is 10, filter out variants observed 11 times or more) in samples from a different cohort/disease target | Yes |
| PLATFORM | Used for setting read group information in the BAM files. e.g. illumina | Yes |
| MEDIAN_COVERAGE_THRESHOLD | The coverage level below which a sample should be failed | Yes |
| TOOLS | Base location of all the tools used by the pipeline |  No |
| SCRIPTS | Directory of the various support scripts that the pipeline uses | No |
| PICARD_HOME | The location of Picard tools | No |
| HTSLIB | The directory in which HTSLIB is installed | No |
| GATK | The directory in which the GenomeAnalysisTK.jar is located | No |
| GATK_LEGACY | Set to 'true' if you are using a version of GATK before version 2.8 | No |
| GATK_VARIANT_ONLY | For versions of GATK before 3.5. If true, do not do genotyping, instead directly call the variants | No |
| EXCEL | Directory in which we store utilities for making Excel files | No |
| BEDTOOLS | bedtools installation directory | No |
| BCFTOOLS | bcftools installation directory | No |
| SAMTOOLS | samtools installation directory | No |
| FASTQC | fastqc installation directory | No |
| GROOVY_NGS | Directory in which the groovy-ngs-utils.jar file is located. Used for processing NGS data in Java/Groovy | No |
| VEP_VERSION | Version of the Ensembl tools from which your VEP installation came, e.g. 83 | No |
| VEP | Directory in which the variant_effect_predictor.pl script is located | No |
| VEP_CACHE | Directory in which the VEP offline cache is located | No |
| IGVTOOLS | Installation directory for IGV tools | No |
| IGV | Installation directory for IGV | No |
| BWA | BWA installation directory | No |
| BWA_THREADS | The number of threads to use when running BWA | Yes |
| CONDEL | Location of the Condel configuration directory | No |
| DBNSFP | DBNSFP reference directory | No |
| VARIANT_DB | Location of the variant database, updated for each sample | No |
| ID_FILE | Location of the ID file, which should contain an ID which uniquely identifies your pipeline installation | No |
| SAMPLE_DB | The central db of all analyses | No |
| UPDATE_VARIANT_DB | By default variant counts are annotated from the same database as the one that they were added to in the first place. This options allows your to specify a separate database for variant annotations | No |
| ANNOTATION_VARIANT_DB | See above | No |
| GROOVY_HOME | Location of Groovy installation | No |
| GROOVY_VERSION | Version of Groovy in GROOVY_HOME | No |
| GROOVY | Location of the Groovy binary | No |
| CHECK_FASTQC_FAILURES | Whether to fail analysis if FASTQC produces warnings | Yes |
| ADAPTERS_FASTA | Set to the adapter sequence for exome capture, if known. Otherwise `false` | Yes |
| EXCLUDE_VARIANT_TYPES | By default synonymous variants are excluded from all outputs. Set to true to include them | Yes |
| JAVA | The path or entry in PATH of the java binary to use | yes
| JAVA_OPTS | Additional flags to use when calling JAVA. e.g. -noverify for certain buggy java builds | Yes |
| SNPEFF | SNPEFF location (not needed by default pipeline) | No |
| HG19_CHROM_INFO | CSV defining the chromosome sizes. Needed for expanded splice regions | No |
| TRIMMOMATIC | Location of the trimmomatic tool (optional) | No |
| BAMSURGEON | Location of the bamsurgeon tool (optional) | No 
| PYTHON | The location or name in PATH of the python executable to use in Cpipe | No |
| splice_region_window | The distance in bp from end of exon from which a variant is annotated as a splicing variant | Yes |
| ANNOTATE_CUSTOM_REGIONS | Bed file to add additional annotation to the final output | Yes |
| INTERVAL_PADDING_CALL | Interval padding to pass to the variant caller | Yes |
| INTERVAL_PADDING_SNV | Padding for SNVs in filter_variants | Yes |
| INTERVAL_PADDING_INDEL | Padding for Indels in filter_variants | Yes |
| ALLOW_SYNONYMOUS_INTRON | Allow synonymous variants this many bases into the intron | Yes |
| ALLOW_SYNONYMOUS_EXON | Allow synonymous variants this many bases into the exon | Yes |
| POST_ANALYSIS_READ_ONLY | Mark the batch read only on completion | Yes |
| POST_ANALYSIS_MOVE | Move the batch directory on completion | Yes
| GAP_ANNOTATOR_CUSTOM_BEDS | Space separated list of additional bed files to consult when generating the gap annotation file | Yes |
| HARD_FILTER_AD | Minimum allele depth | Yes |
| HARD_FILTER_AF | Minimum allele frequency | Yes |
| HARD_FILTER_DP | Minimum depth | Yes |
| HARD_FILTER_QUAL | Minimum quality | Yes |
| TRIO_REFINEMENT_SUPPORTING | Location of the trio refinement file | No |
| READ_PERCENTAGE_THRESHOLD | |
| FILTERED_ON_EXONS | Generate a filtered bam file. Either exons, design, or skip | Yes |
| QC_THRESHOLD | What depth is required to contribute to satisfactory coverage | Yes |
| QC_GOOD | What percentage of QC_THRESHOLD must be achieved across the gene to get a good rating | Yes |
| QC_PASS | What percentage of QC_THRESHOLD must be achieved across the gene to get a pass rating | Yes |
| QC_FAIL | What percentage of QC_THRESHOLD must be achieved across the gene to get a fail rating | Yes |

## bpipe.config
The bpipe configuration file is a single file named `bpipe.config` located in the pipeline directory, which controls 
some aspects of bpipe's functionality. You can find out more information about this file from reading the bpipe docs
 [available here](https://github.com/ssadedin/bpipe/tree/master/docs). Specifically, `bpipe.config` controls:
 * Running of bpipe with Cluster Management Systems (TORQUE, Slurm etc.)
 * [Sending of status emails for the pipeline](https://github.com/ssadedin/bpipe/blob/master/docs/Guides/Notifications.md#notifications-in-bpipe)
 * [Specifying the amount of concurrency with which the pipeline is run](https://github.com/ssadedin/bpipe/blob/master/docs/Language/ParallelTasks.md)
 * [Concurrency for scanning output files](https://github.com/ssadedin/bpipe/blob/master/docs/Guides/AdvancedSettings.md#output-scan-concurrency)
 * [Separation time between jobs](https://github.com/ssadedin/bpipe/blob/master/docs/Guides/AdvancedSettings.md#job-launch-separation) (stages of the pipeline)
 * [The location of a custom R executable to be used by bpipe](https://github.com/ssadedin/bpipe/blob/master/docs/Language/R.md#executing-inline-r-code) (cpipe has its own internal version)
