# Terminology

On this page we explain some of the terms we use in this documentation to describe the pipline and analysis

* Sample: A single fastq sequence or pair of fastqs that come from the same patient, intended to be analysed
* Design: Also called a cohort or profile. A set of genes or regions that are relevant to a particular condition. Samples
analysed with a given design only consider those regions and ignore variance in all other regions
* Batch: A group of samples designed to be analysed at the same time
* Profile: See 'Design'
* Cohort: See 'Design'
* Sample Metadata: Data about the samples, e.g. sample ID, sex, ethnicity
* Pipeline Run ID: A number, designed to be unique to each run of cpipe across any installations on any system. Generated
using the string contained in the `pipeline_id` file.