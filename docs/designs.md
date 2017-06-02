# Designs

* [Introduction](#introduction)
* [Gene List](#gene-list)
* [Regions BED file](#regions-bed-file)

## Introduction

In Cpipe, a design specifies the analysed regions of the genome for a given set of samples. For instance, a given disease might
associated with a group of genes, which we then group into a design and re-use for analysis of patients being tested
for this disease.
 
However, our current recommendation is that you use the existing ALL design, which contains all genes in the UCSC database
as this will save you having to re-analyse the data in case your gene list changes.
 
The design used to analyse a given sample is specified in the cohort column of the sample metadata file. See the [sample 
metadata documentation](batches.md#sample-metadata) for more information.

Cpipe also contains a number of useful command line tools for editing designs. Refer to the documentation 
[here](commands.md#design)
 
 
A design consists of the following files
* [A gene list (mandatory)](#gene-list)
* [A regions BED file (optional)](#regions-bed-file)

## Gene List
A gene list is expressed as a tab-separated (TSV) text file named `<profile name>.genes.txt`, located in the design directory 
(designs/<profile name>). The gene list has the following fields per line
* The offical HGNC gene symbol 
* The priority of this gene, from 1-4, where 1 is the lowest and 4 is the highest priority

If the gene list is specified and not a BED file (see below), the gene list is expanded out to a BED file describing the 
regions encompassing the genes of interest using the UCSC refseq exon database. 
This conversion is performed automatically by Cpipe.
 
## Regions BED file
If a BED named `<profile name>.bed` is present in the profile root (`designs/<profile name>`), then it will be used
to determine the analysed regions instead of the gene list. However, the gene list must still be present in order to
generate QC reports, so this cannot be used entirely instead of a gene list.