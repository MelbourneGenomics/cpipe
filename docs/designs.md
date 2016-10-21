# Designs

## Introduction

In Cpipe, a design specifies the analysed regions of the genome for a given set of samples. For instance, a given disease might
 associated with a group of genes, which we then group into a design and re-use for analysis of patients being tested
 for this disease.
 
A design consists of the following files
* [A gene list (mandatory)](#gene-list)
* [A regions BED file (optional)](#regions-bed-file)

## Gene List
A gene list is expressed as a tab-separated (TSV) text file named '<profile name>.genes.txt', located in the design directory 
(designs/<profile name>). The gene list has the following fields per line
* The offical HGNC gene symbol 
* The priority of this gene, from 1-4, where 1 is the lowest and 4 is the highest priority

If the gene list is specified and not a BED file (see below), the gene list is expanded out to a BED file describing the 
regions encompassing the genes of interest using the UCSC refseq exon database. 
This conversion is performed automatically by Cpipe.
 
## Regions BED file
If a BED named <profile name>.bed is present in the profile root (designs/<profile name>), then it will be used
to determine the analysed regions instead of the gene list. However, the gene list must still be present in order to
generate QC reports.