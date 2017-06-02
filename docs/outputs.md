# Outputs

* [Batch Validation Files](#batch-validation-files)
  * [Gender Validation](#gender-validation)
  * [Gene Coverage By Sample](#gene-coverage-by-sample)
  * [Observed mean coverage by sample](#observed-mean-coverage-by-sample)
  * [Individual genes with &gt;80\x fail across samples](#individual-genes-with-80-fail-across-samples)
  * [Requested Genes not found in Reference](#requested-genes-not-found-in-reference)
  * [Requested Genes not found in Annovar](#requested-genes-not-found-in-annovar)
* [Summary Report](#summary-report)
  * [Summary Data](#summary-data)
    * [Inferred Sex Calculation](#inferred-sex-calculation)
  * [Coverage Summary](#coverage-summary)
    * [Mean Coverage Calculation](#mean-coverage-calculation)
  * [Gene Summary](#gene-summary)
* [Tab Separated Variant File](#tab-separated-variant-file)
* [Provenance Report](#provenance-report)
* [Gap Report](#gap-report)

***

## Introduction

After completion of the pipeline (after a successful `run` command), cpipe will produce results in the `cpipe/batches/
<batch identifier>/analysis` directory. This directory has the following important subdirectories:

* `align`: Contains the sequence alignment (BAM) files generated when aligning the fastq sample files to the 
reference genome. 
* `fastqc`: Contains quality control information on the fastq sample files
* `qc`: Contains miscellaneous quality control on the sequence alignment
* `results`: The main results directory. Contains all of the final cpipe outputs. The important files in this directory
 are these (documented below):
    * [Batch validation files](#batch-validation-files)
    * [Summmary reports](#summary-report)
    * [Tab Separated Variant File](#tab-separated-variant-file)
    * [Provenance Report](#provenance-report)
    * [Gap Report](#gap-report)
* `variants`: Contains the raw variant (VCF) files. The VCF files that are likely to be the most interesting are:
    * `<sample ID>.individual.genotype.soi.vep.post_filter.vcfanno.vcf`: The last and most completely annotated and
    normalised VCF. Use this VCF for importing into curation tools etc.
    * `<sample ID>.individual.refined.vcf`: The raw VCF that has been run through GATK but has not been normalised or 
    annotated.

Note that the following types of outputs have variable file names, all of which include the pipeline run ID. If you wish
to understand the meaning of this, refer to the [terminology](./terminology.md) section.

## Batch Validation Files
The batch validation report summarizes common sanity checks performed on a batch to quickly verify that sample results
appear valid.

These files are named `<pipeline run ID>_batch_validation.html/md`. The .md markdown files are designed to be read as 
plain text (e.g. in a terminal), while the .html files are designed to be viewed graphically (e.g. from a web browser),
although they present the same content.

The report consists of the following sections:

### Gender Validation
This result is taken from the sample’s `summary.karyotype.tsv` file.

### Gene Coverage By Sample
* This result is a summary of the sample’s summary.pdf file, which calculates if a gene has acceptable coverage.
* The number of genes that are good, OK, and failed is counted and displayed.
* If the proportion of failed genes is greater than 15% then the outcome for the sample is marked as fail.

### Observed mean coverage by sample
* This result is based on the observed mean coverage displayed in the sample’s summary.pdf file.
* The result for each sample is displayed, with an observed mean coverage of less than 90 being marked as fail.

### Individual genes with >80% fail across samples
We look at coverage results by gene rather than sample, aggregating all outcomes for a gene across the batch. 
If a gene has failed to be acceptably covered in >80% of cases, it is displayed here as a fail.

### Requested Genes not found in Reference
Genes that were requested in the cohort or in the sample metadata file that are not listed in the refseq genes used to generate the bed file.

### Requested Genes not found in Annovar
Genes that were requested in the cohort or in the sample metadata file that are not found in Annovar’s list of genes 
(deprecated since Annovar is no longer used by Cpipe).

## Summary Report
Cpipe generates a summary HTML for each analysed sample in a batch. It is of the form 
`<sample_run_id>_<sample_name>.summary.htm`(or `.md`).
 
The report contains three main sections:
* [Summary Data](#summary-data)
* [Coverage Summary](#coverage-summary)
* [Gene Summary](#gene-summary)

### Summary Data

This section contains the following fields:
 
| Field Name | Source |
| --- | --- |
| Batch | sample metadata file |
| Study ID | sample metadata file |
| Sex | sample metadata file |
| Inferred Sex | [See detail below](#inferred-sex-calculation).  Displayed in red if this value doesn’t match Sex.|
| Disease Cohort | sample metadata file |
| Hospital/Institution | sample metadata file |
| Ethnicity | sample metadata file |
| Prioritized Genes | sample metadata file |
| Consanguinity Status | sample metadata file (with _ replaced with spaces) |
| Sample type (tumor/normal) | sample metadata file | 
| Sequencing Dates | sample metadata file |
| DNA Collection Dates | sample metadata file |
| Sequencing Machines | sample metadata file |
 
#### Inferred Sex Calculation

We use the coverage from chromosomes 1 and 22 as proxy estimates for the autosomal coverage.  
For females we expect the X coverage to be similar to the autosomes coverage (i.e. diploid).  
For males, we expect the X coverage to be less than the autosome coverage.

* The mean coverage on chromosome 1 and chromosome 22 is calculated.
* The mean coverage on chromosome X is calculated.
* The mean coverage on chromosome Y is calculated.

For females we expect the chr Y coverage to be low, i.e. less than 5X, and the chr X coverage to be > 30X.
These figures assume Nextera exome capture with at least 100X coverage overall and this Female inference would need to
be expanded if significantly higher or lower overall coverage is expected from the WES test.
To ensure robustness to overall test coverage, we might consider developing a ratio based metric for inferring Female 
sex (i.e. chr X mean coverage / Chr Y mean coverage > 3 )

* If chrY mean coverage < 5 and chrX mean coverage > 30 then the inferred sex is FEMALE
* If the chr Y mean coverage is > 5X, this could indicate a Male, or it could be a Female with very high overall coverage (such that the chr Y reaches above 5X by chance).  Therefore, to check for a Male, we base the test on the relative mean coverage of chr X with the autosomes.  We expect a ratio of a half for Males (hence the cut-off of 0.7), while Females we expect 1.
* Otherwise, if mean coverage across chrX / mean coverage across chr1, chr22 < 0.7 then the inferred sex is MALE
* Otherwise, the inferred sex is OTHER

We note that in Nextera captures, the chr Y coverage seems to be higher than expected, often achieving similar mean coverage to the autosomes in Males.  It is for this reason that the sex inference for Males does not use the ratio of chrY coverage compared to the autosome.

### Coverage Summary

The coverage summary includes the following fields:

* Mean Coverage Reported by Lab: the mean coverage reported by the sequencing lab
* Observed Mean Coverage: the mean coverage across the capture region
* Observed Median Coverage: the median coverage across the capture region
* Total Reads: the total number of reads generated by the sequencer
* Unmapped Reads: reads that were not mapped to the genome
* Mapped Paired Reads: paired reads that were mapped to the genome
* % Mapped on Target: the proportion of mapped reads that have any part align to any part of the capture region
* % Coverage within 20% of Mean: bases in the capture region with coverage within 20% of the observed mean coverage
* Mean Fragment Size: the average distance between correctly mapped and paired reads
* Perc: the percentage of the gene overlapping the capture region with acceptable coverage
* Median: the median coverage across the gene overlapping the capture region
* % in capture: the proportion of the gene that overlaps the capture region

#### Mean Coverage Calculation
* Regions covered by the exome capture region are considered.
* Only mapped reads with a mapping quality of 1 or more are considered.
* All mapped reads that overlap a region are considered
* Mean Coverage = number of reads overlapping each position in each region / total number of positions

### Gene Summary
This section displays data for a list of genes. The list of genes comes from the coverage file, which ultimately comes 
from the disease cohort BED file.

The table has the following fields:

| Field Name | Source |
| --- | --- |
| Gene | The gene name from the coverage file |
| Category | category from cohort.genes.txt |
| Perc > 20x | Percentage of the bases of this gene (in the cohort BED) with coverage greater than 20 |
| Median | The median coverage over this gene (See below definition for gene coverage) |
| OK? | Based on “Perc > 20X”: >95%: Good, >80%: Pass, Otherwise: Fail |

## Tab Separated Variant File
A compilation of all the variant-level data is found in the files named 
`<pipeline run ID>_<sample ID>.<individual|trio>.lovd.tsv`. This file has one row per 
variant in the individual, with the following columns:

ID | Description | Present in Individual | Present in Trio | Imported into LOVD | Annotation Source
--- | --- | --- | --- | --- | --- 
AD | Allelic depths for the ref and alt alleles in the order listed | Yes | Yes | Yes | GATK
DP | Approximate read depth (reads with MQ=255 or with bad mates are filtered) | Yes | Yes | Yes | GATK
FT | Genotype-level filter | No | Yes | No | GATK
GQ | Genotype Quality | Yes | Yes | Yes | GATK
GT | Genotype | Yes | Yes | Yes | GATK
JL | Phred-scaled joint likelihood of the genotype combination (before applying family priors)  | Yes | No | Yes | GATK
JP | Phred-scaled joint posterior probability of the genotype combination (after applying family priors)  | Yes | No | Yes | GATK
MIN_DP | Minimum DP observed within the GVCF block | Yes | Yes | No | GATK
PGT | Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another | Yes | Yes | Yes | GATK
PID | Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group | Yes | Yes | Yes | GATK
PL | Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification | Yes | Yes | Yes | GATK
PP | Phred-scaled Posterior Genotype Probabilities | No | Yes | Yes | GATK
RGQ | Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong) | Yes | Yes | No | GATK
SB | Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias. | Yes | Yes | No | GATK
ABHet | Allele Balance for heterozygous calls (ref/(ref+alt)) | Yes | Yes | Yes | GATK
ABHom | Allele Balance for homozygous calls (A/(A+O)) where A is the allele (ref or alt) and O is anything other | Yes | Yes | Yes | GATK
AC | Allele count in genotypes, for each ALT allele, in the same order as listed | Yes | Yes | Yes | GATK
AF | Allele Frequency, for each ALT allele, in the same order as listed | Yes | Yes | Yes | GATK
AN | Total number of alleles in called genotypes | Yes | Yes | Yes | GATK
BaseQRankSum | Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities | Yes | Yes | Yes | GATK
ClippingRankSum | Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases | Yes | Yes | Yes | GATK
DB | dbSNP Membership | Yes | Yes | Yes | GATK
DP | Approximate read depth; some reads may have been filtered | Yes | Yes | Yes | GATK
DS | Were any of the samples downsampled? | Yes | Yes | No | GATK
END | Stop position of the interval | Yes | Yes | No | GATK
ExcessHet | Phred-scaled p-value for exact test of excess heterozygosity | Yes | Yes | Yes | GATK
FS | Phred-scaled p-value using Fisher's exact test to detect strand bias | Yes | Yes | Yes | GATK
GC | GC content around the variant (see docs for window size details) | Yes | Yes | No | GATK
GQ_MEAN | Mean of all GQ values | Yes | Yes | Yes | GATK
GQ_STDDEV | Standard deviation of all GQ values | Yes | Yes | No | GATK
HaplotypeScore | Consistency of the site with at most two segregating haplotypes | Yes | Yes | No | GATK
InbreedingCoeff | Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation | Yes | Yes | No | GATK
LikelihoodRankSum | Z-score from Wilcoxon rank sum test of Alt Vs. Ref haplotype likelihoods | Yes | Yes | Yes | GATK
MLEAC | Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed | Yes | Yes | Yes | GATK
MLEAF | Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed | Yes | Yes | Yes | GATK
MQ | RMS Mapping Quality | Yes | Yes | Yes | GATK
MQRankSum | Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities | Yes | Yes | Yes | GATK
NCC | Number of no-called samples | Yes | Yes | No | GATK
OND | Overall non-diploid ratio (alleles/(alleles+non-alleles)) | Yes | Yes | Yes | GATK
PG | Genotype Likelihood Prior  | No | Yes | No | GATK
QD | Variant Confidence/Quality by Depth | Yes | Yes | Yes | GATK
RAW_MQ | Raw data for RMS Mapping Quality | Yes | Yes | No | GATK
ReadPosRankSum | Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias | Yes | Yes | Yes | GATK
SOR | Symmetric Odds Ratio of 2x2 contingency table to detect strand bias | Yes | Yes | Yes | GATK
VariantType | Variant type description | Yes | Yes | Yes | GATK
hiConfDeNovo | High confidence possible de novo mutation (GQ >= 20 for all trio members)=[comma-delimited list of child samples]  | No | Yes | Yes | GATK
loConfDeNovo | Low confidence possible de novo mutation (GQ >= 10 for child, GQ > 0 for parents)=[comma-delimited list of child samples]  | No | Yes | Yes | GATK
set | Source VCF for the merged record in CombineVariants | Yes | Yes | Yes | GATK
Allele | Variant allele used to calculate the consequence | Yes | Yes | Yes | VEP
Consequence | Consequence type of this variant | Yes | Yes | Yes | VEP
IMPACT | Impact modifier for the consequence type | Yes | Yes | Yes | VEP
SYMBOL | Gene symbol | Yes | Yes | Yes | VEP
Gene | Ensembl stable ID of affected gene | Yes | Yes | Yes | VEP
Feature_type | Type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature. | Yes | Yes | Yes | VEP
Feature | Ensembl stable ID of feature | Yes | Yes | Yes | VEP
BIOTYPE | Biotype of transcript or regulatory feature | Yes | Yes | Yes | VEP
EXON | Exon number (out of total number) | Yes | Yes | Yes | VEP
INTRON | Intron number (out of total number) | Yes | Yes | Yes | VEP
HGVSc | HGVS coding sequence name | Yes | Yes | Yes | VEP
HGVSp | HGVS protein sequence name | Yes | Yes | Yes | VEP
cDNA_position | Relative position of base pair in cDNA sequence | Yes | Yes | Yes | VEP
CDS_position | Relative position of base pair in coding sequence | Yes | Yes | Yes | VEP
Protein_position | Relative position of amino acid in protein | Yes | Yes | Yes | VEP
Amino_acids | Only given if the variant affects the protein-coding sequence | Yes | Yes | Yes | VEP
Codons | The alternative codons with the variant base in upper case | Yes | Yes | Yes | VEP
Existing_variation | Known variants | Yes | Yes | Yes | VEP
ALLELE_NUM | Allele number from input; 0 is reference, 1 is first alternate etc | Yes | Yes | Yes | VEP
DISTANCE | Shortest distance from variant to transcript | Yes | Yes | Yes | VEP
STRAND | DNA strand (1 or -1) on which the transcript/feature lies | Yes | Yes | Yes | VEP
SYMBOL_SOURCE | Source of the gene symbol | Yes | Yes | Yes | VEP
HGNC_ID | HGNC ID | Yes | Yes | Yes | VEP
CANONICAL | A flag indicating if the transcript is denoted as the canonical transcript for this gene | Yes | Yes | Yes | VEP
ENSP | Ensembl protein identifier of the affected transcript | Yes | Yes | Yes | VEP
REFSEQ_MATCH | RefSeq transcript match status; contains a number of flags indicating whether this RefSeq transcript matches the underlying reference sequence and/or an Ensembl transcript: | Yes | Yes | Yes | VEP
SIFT | SIFT prediction and/or score, with both given as prediction(score) | Yes | Yes | Yes | VEP
PolyPhen | PolyPhen prediction and/or score | Yes | Yes | Yes | VEP
HGVS_OFFSET | Indicates by how many bases the HGVS notations for this variant have been shifted | Yes | Yes | Yes | VEP
GMAF | Non-reference allele and frequency of existing variant in 1000 Genomes | Yes | Yes | Yes | VEP
AFR_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined African population | Yes | Yes | Yes | VEP
AMR_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined American population | Yes | Yes | Yes | VEP
EAS_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined East Asian population | Yes | Yes | Yes | VEP
EUR_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined European population | Yes | Yes | Yes | VEP
SAS_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined South Asian population | Yes | Yes | Yes | VEP
AA_MAF | Non-reference allele and frequency of existing variant in NHLBI-ESP African American population | Yes | Yes | Yes | VEP
EA_MAF | Non-reference allele and frequency of existing variant in NHLBI-ESP European American population | Yes | Yes | Yes | VEP
ExAC_MAF | Frequency of existing variant in ExAC combined population | Yes | Yes | Yes | VEP
ExAC_Adj_MAF | Adjusted frequency of existing variant in ExAC combined population | Yes | Yes | Yes | VEP
ExAC_AFR_MAF | Frequency of existing variant in ExAC African/American population | Yes | Yes | Yes | VEP
ExAC_AMR_MAF | Frequency of existing variant in ExAC American population | Yes | Yes | Yes | VEP
ExAC_EAS_MAF | Frequency of existing variant in ExAC East Asian population | Yes | Yes | Yes | VEP
ExAC_FIN_MAF | Frequency of existing variant in ExAC Finnish population | Yes | Yes | Yes | VEP
ExAC_NFE_MAF | Frequency of existing variant in ExAC Non-Finnish European population | Yes | Yes | Yes | VEP
ExAC_OTH_MAF | Frequency of existing variant in ExAC combined other combined populations | Yes | Yes | Yes | VEP
ExAC_SAS_MAF | Frequency of existing variant in ExAC South Asian population | Yes | Yes | Yes | VEP
CLIN_SIG | Clinical significance of variant from dbSNP | Yes | Yes | Yes | VEP
SOMATIC | Somatic status of existing variant(s) | Yes | Yes | Yes | VEP
PHENO | Indicates if existing variant is associated with a phenotype, disease or trait | Yes | Yes | Yes | VEP
PUBMED | Pubmed ID(s) of publications that cite existing variant | Yes | Yes | Yes | VEP
1000Gp1_AC | Alternative allele counts in the whole 1000 genomes phase 1 (1000Gp1) data. | Yes | Yes | Yes | DbNSFP VEP Plugin
1000Gp1_AF | Alternative allele frequency in the whole 1000Gp1 data. | Yes | Yes | Yes | DbNSFP VEP Plugin
1000Gp1_AFR_AC | Alternative allele counts in the 1000Gp1 African descendent samples. | Yes | Yes | Yes | DbNSFP VEP Plugin
1000Gp1_AFR_AF | Alternative allele frequency in the 1000Gp1 African descendent samples. | Yes | Yes | Yes | DbNSFP VEP Plugin
1000Gp1_AMR_AC | Alternative allele counts in the 1000Gp1 American descendent samples. | Yes | Yes | No | DbNSFP VEP Plugin
1000Gp1_AMR_AF | Alternative allele frequency in the 1000Gp1 American descendent samples. | Yes | Yes | No | DbNSFP VEP Plugin
1000Gp1_ASN_AC | Alternative allele counts in the 1000Gp1 Asian descendent samples. | Yes | Yes | Yes | DbNSFP VEP Plugin
1000Gp1_ASN_AF | Alternative allele frequency in the 1000Gp1 Asian descendent samples. | Yes | Yes | No | DbNSFP VEP Plugin
1000Gp1_EUR_AC | Alternative allele counts in the 1000Gp1 European descendent samples. | Yes | Yes | No | DbNSFP VEP Plugin
1000Gp1_EUR_AF | Alternative allele frequency in the 1000Gp1 European descendent samples. | Yes | Yes | Yes | DbNSFP VEP Plugin
CADD_phred | CADD phred-like score. This is phred-like rank score based on whole genome CADD raw scores. Please refer to Kircher et al. (2014) Nature Genetics 46(3):310-5 for details. The larger the score the more likely the SNP has damaging effect. Please note the following copyright statement for CADD: "CADD scores (http://cadd.gs.washington.edu/) are Copyright 2013 University of Washington and Hudson-Alpha Institute for Biotechnology (all rights reserved) but are freely available for all academic, non-commercial applications. For commercial licensing information contact Jennifer McCullar (mccullaj@uw.edu)." | Yes | Yes | Yes | DbNSFP VEP Plugin
CADD_raw | CADD raw score for funtional prediction of a SNP. Please refer to Kircher et al. (2014) Nature Genetics 46(3):310-5 for details. The larger the score the more likely the SNP has damaging effect. Please note the following copyright statement for CADD: "CADD scores (http://cadd.gs.washington.edu/) are Copyright 2013 University of Washington and Hudson-Alpha Institute for Biotechnology (all rights reserved) but are freely available for all academic, non-commercial applications. For commercial licensing information contact Jennifer McCullar (mccullaj@uw.edu)." | Yes | Yes | No | DbNSFP VEP Plugin
CADD_raw_rankscore | CADD raw scores were ranked among all CADD raw scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of CADD raw scores in dbNSFP. Please note the following copyright statement for CADD: "CADD scores (http://cadd.gs.washington.edu/) are Copyright 2013 University of Washington and Hudson-Alpha Institute for Biotechnology (all rights reserved) but are freely available for all academic, non-commercial applications. For commercial licensing information contact Jennifer McCullar (mccullaj@uw.edu)." | Yes | Yes | Yes | DbNSFP VEP Plugin
ESP6500_AA_AF | Alternative allele frequency in the Afrian American samples of the NHLBI GO Exome Sequencing Project (ESP6500 data set). | Yes | Yes | Yes | DbNSFP VEP Plugin
ESP6500_EA_AF | Alternative allele frequency in the European American samples of the NHLBI GO Exome Sequencing Project (ESP6500 data set). | Yes | Yes | Yes | DbNSFP VEP Plugin
FATHMM_pred | If a FATHMMori score is <=-1.5 (or rankscore <=0.81415) the corresponding NS is predicted as "D(AMAGING)"; otherwise it is predicted as "T(OLERATED)". Multiple predictions separated by ";" | Yes | Yes | Yes | DbNSFP VEP Plugin
FATHMM_rankscore | FATHMMori scores were ranked among all FATHMMori scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of FATHMMori scores in dbNSFP. If there are multiple scores, only the most damaging (largest) rankscore is presented. The scores range from 0 to 1. | Yes | Yes | Yes | DbNSFP VEP Plugin
FATHMM_score | FATHMM default score (weighted for human inherited-disease mutations with Disease Ontology) (FATHMMori). Scores range from -18.09 to 11.0. Multiple scores separated by ";" Please refer to Shihab et al. (2013) Human Mutation 34(1):57-65 for details. | Yes | Yes | Yes | DbNSFP VEP Plugin
GERP++_NR | GERP++ neutral rate | Yes | Yes | Yes | DbNSFP VEP Plugin
GERP++_RS | GERP++ RS score, the larger the score, the more conserved the site. | Yes | Yes | Yes | DbNSFP VEP Plugin
GERP++_RS_rankscore | GERP++ RS scores were ranked among all GERP++ RS scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of GERP++ RS scores in dbNSFP. | Yes | Yes | Yes | DbNSFP VEP Plugin
LRT_Omega | estimated nonsynonymous-to-synonymous-rate ratio (Omega, reported by LRT) | Yes | Yes | Yes | DbNSFP VEP Plugin
LRT_converted_rankscore | LRTori scores were first converted as LRTnew=1-LRTori*0.5 if Omega<1, or LRTnew=LRTori*0.5 if Omega>=1. Then LRTnew scores were ranked among all LRTnew scores in dbNSFP. The rankscore is the ratio of the rank over the total number of the scores in dbNSFP. The scores range from 0.00166 to 0.85682. | Yes | Yes | Yes | DbNSFP VEP Plugin
LRT_pred | LRT prediction, D(eleterious), N(eutral) or U(nknown), which is not solely determined by the score. | Yes | Yes | Yes | DbNSFP VEP Plugin
LRT_score | The original LRT two-sided p-value (LRTori), ranges from 0 to 1. | Yes | Yes | Yes | DbNSFP VEP Plugin
MetaLR_pred | Prediction of our MetaLR based ensemble prediction score,"T(olerated)" or "D(amaging)". The score cutoff between "D" and "T" is 0.5. The rankscore cutoff between "D" and "T" is 0.82268. | Yes | Yes | Yes | DbNSFP VEP Plugin
MetaLR_rankscore | MetaLR scores were ranked among all MetaLR scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of MetaLR scores in dbNSFP. The scores range from 0 to 1. | Yes | Yes | Yes | DbNSFP VEP Plugin
MetaLR_score | Our logistic regression (LR) based ensemble prediction score, which incorporated 10 scores (SIFT, PolyPhen-2 HDIV, PolyPhen-2 HVAR, GERP++, MutationTaster, Mutation Assessor, FATHMM, LRT, SiPhy, PhyloP) and the maximum frequency observed in the 1000 genomes populations. Larger value means the SNV is more likely to be damaging. Scores range from 0 to 1. | Yes | Yes | Yes | DbNSFP VEP Plugin
MetaSVM_pred | Prediction of our SVM based ensemble prediction score,"T(olerated)" or "D(amaging)". The score cutoff between "D" and "T" is 0. The rankscore cutoff between "D" and "T" is 0.83357. | Yes | Yes | Yes | DbNSFP VEP Plugin
MetaSVM_rankscore | MetaSVM scores were ranked among all MetaSVM scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of MetaSVM scores in dbNSFP. The scores range from 0 to 1. | Yes | Yes | Yes | DbNSFP VEP Plugin
MetaSVM_score | Our support vector machine (SVM) based ensemble prediction score, which incorporated 10 scores (SIFT, PolyPhen-2 HDIV, PolyPhen-2 HVAR, GERP++, MutationTaster, Mutation Assessor, FATHMM, LRT, SiPhy, PhyloP) and the maximum frequency observed in the 1000 genomes populations. Larger value means the SNV is more likely to be damaging. Scores range from -2 to 3 in dbNSFP. | Yes | Yes | Yes | DbNSFP VEP Plugin
MutationAssessor_pred | MutationAssessor's functional impact of a variant : predicted functional, i.e. high ("H") or medium ("M"), or predicted non-functional, i.e. low ("L") or neutral ("N"). The MAori score cutoffs between "H" and "M", "M" and "L", and "L" and "N", are 3.5, 1.935 and 0.8, respectively. The rankscore cutoffs between "H" and "M", "M" and "L", and "L" and "N", are 0.92924, 0.51945 and 0.19692, respectively. | Yes | Yes | Yes | DbNSFP VEP Plugin
MutationAssessor_rankscore | MAori scores were ranked among all MAori scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of MAori scores in dbNSFP. The scores range from 0 to 1. | Yes | Yes | Yes | DbNSFP VEP Plugin
MutationAssessor_score | MutationAssessor functional impact combined score (MAori). The score ranges from -5.135 to 6.49 in dbNSFP. Please refer to Reva et al. (2011) Nucl. Acids Res. 39(17):e118 for details. | Yes | Yes | Yes | DbNSFP VEP Plugin
MutationTaster_converted_rankscore | The MTori scores were first converted: if the prediction is "A" or "D" MTnew=MTori; if the prediction is "N" or "P", MTnew=1-MTori. Then MTnew scores were ranked among all MTnew scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of MTnew scores in dbNSFP. The scores range from 0.09067 to 0.80722. | Yes | Yes | Yes | DbNSFP VEP Plugin
MutationTaster_pred | MutationTaster prediction, "A" ("disease_causing_automatic"), "D" ("disease_causing"), "N" ("polymorphism") or "P" ("polymorphism_automatic"). The score cutoff between "D" and "N" is 0.5 for MTori and 0.31655 for the rankscore. | Yes | Yes | Yes | DbNSFP VEP Plugin
MutationTaster_score | MutationTaster p-value (MTori), ranges from 0 to 1. | Yes | Yes | Yes | DbNSFP VEP Plugin
PROVEAN_converted_rankscore | PROVEANori were first converted to PROVEANnew=1-(PROVEANori+14)/28, then ranked among all PROVEANnew scores in dbNSFP. The rankscore is the ratio of the rank the PROVEANnew score over the total number of PROVEANnew scores in dbNSFP. If there are multiple scores, only the most damaging (largest) rankscore is presented. | Yes | Yes | Yes | DbNSFP VEP Plugin
PROVEAN_pred | If PROVEANori <= -2.5 (rankscore>=0.59) the corresponding NS is predicted as "D(amaging)"; otherwise it is predicted as "N(eutral)". Multiple predictions separated by ";" | Yes | Yes | Yes | DbNSFP VEP Plugin
PROVEAN_score | PROVEAN score (PROVEANori). Scores range from -14 to 14. The smaller the score the more likely the SNP has damaging effect. Multiple scores separated by ";". Details can be found in DOI: 10.1371/journal.pone.0046688 | Yes | Yes | Yes | DbNSFP VEP Plugin
Polyphen2_HDIV_pred | Polyphen2 prediction based on HumDiv, "D" ("porobably damaging", HDIV score in [0.957,1] or rankscore in [0.52996,0.89917]), "P" ("possibly damaging", HDIV score in [0.453,0.956] or rankscore in [0.34412,0.52842]) and "B" ("benign", HDIV score in [0,0.452] or rankscore in [0.02656,0.34399]). Score cutoff for binary classification is 0.5 for HDIV score or 0.35411 for rankscore, i.e. the prediction is "neutral" if the HDIV score is smaller than 0.5 (rankscore is smaller than 0.35411), and "deleterious" if the HDIV score is larger than 0.5 (rankscore is larger than 0.35411). Multiple entries are separated by ";". | Yes | Yes | Yes | DbNSFP VEP Plugin
Polyphen2_HDIV_rankscore | Polyphen2 HDIV scores were first ranked among all HDIV scores in dbNSFP. The rankscore is the ratio of the rank the score over the total number of the scores in dbNSFP. If there are multiple scores, only the most damaging (largest) rankscore is presented. The scores range from 0.02656 to 0.89917. | Yes | Yes | Yes | DbNSFP VEP Plugin
Polyphen2_HDIV_score | Polyphen2 score based on HumDiv, i.e. hdiv_prob. The score ranges from 0 to 1. Multiple entries separated by ";". | Yes | Yes | Yes | DbNSFP VEP Plugin
Polyphen2_HVAR_pred | Polyphen2 prediction based on HumVar, "D" ("porobably damaging", HVAR score in [0.909,1] or rankscore in [0.62955,0.9711]), "P" ("possibly damaging", HVAR in [0.447,0.908] or rankscore in [0.44359,0.62885]) and "B" ("benign", HVAR score in [0,0.446] or rankscore in [0.01281,0.44315]). Score cutoff for binary classification is 0.5 for HVAR score or 0.45998 for rankscore, i.e. the prediction is "neutral" if the HVAR score is smaller than 0.5 (rankscore is smaller than 0.45998), and "deleterious" if the HVAR score is larger than 0.5 (rankscore is larger than 0.45998). Multiple entries are separated by ";". | Yes | Yes | Yes | DbNSFP VEP Plugin
Polyphen2_HVAR_rankscore | Polyphen2 HVAR scores were first ranked among all HVAR scores in dbNSFP. The rankscore is the ratio of the rank the score over the total number of the scores in dbNSFP. If there are multiple scores, only the most damaging (largest) rankscore is presented. The scores range from 0.01281 to 0.9711. | Yes | Yes | Yes | DbNSFP VEP Plugin
Polyphen2_HVAR_score | Polyphen2 score based on HumVar, i.e. hvar_prob. The score ranges from 0 to 1. Multiple entries separated by ";". | Yes | Yes | Yes | DbNSFP VEP Plugin
Reliability_index | Number of observed component scores (except the maximum frequency in the 1000 genomes populations) for MetaSVM and MetaLR. Ranges from 1 to 10. As MetaSVM and MetaLR scores are calculated based on imputed data, the less missing component scores, the higher the reliability of the scores and predictions. | Yes | Yes | Yes | DbNSFP VEP Plugin
SIFT_pred | If SIFTori is smaller than 0.05 (rankscore>0.55) the corresponding NS is predicted as "D(amaging)"; otherwise it is predicted as "T(olerated)". Multiple predictions separated by ";" | Yes | Yes | Yes | DbNSFP VEP Plugin
SiPhy_29way_logOdds | SiPhy score based on 29 mammals genomes. The larger the score, the more conserved the site. | Yes | Yes | Yes | DbNSFP VEP Plugin
SiPhy_29way_logOdds_rankscore | SiPhy_29way_logOdds scores were ranked among all SiPhy_29way_logOdds scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of SiPhy_29way_logOdds scores in dbNSFP. | Yes | Yes | Yes | DbNSFP VEP Plugin
SiPhy_29way_pi | The estimated stationary distribution of A, C, G and T at the site, using SiPhy algorithm based on 29 mammals genomes. | Yes | Yes | Yes | DbNSFP VEP Plugin
UniSNP_ids | rs numbers from UniSNP, which is a cleaned version of dbSNP build 129, in format: rs number1;rs number2;... | Yes | Yes | Yes | DbNSFP VEP Plugin
VEST3_rankscore | VEST3 scores were ranked among all VEST3 scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of VEST3 scores in dbNSFP. The scores range from 0 to 1. Please note VEST score is free for non-commercial use. For more details please refer to http://wiki.chasmsoftware.org/index.php/SoftwareLicense. Commercial users should contact the Johns Hopkins Technology Transfer office. | Yes | Yes | Yes | DbNSFP VEP Plugin
VEST3_score | VEST 3.0 score. Score ranges from 0 to 1. The larger the score the more likely the mutation may cause functional change. In case there are multiple scores for the same variant, the largest score (most damaging) is presented. Please refer to Carter et al., (2013) BMC Genomics. 14(3) 1-16 for details. Please note this score is free for non-commercial use. For more details please refer to http://wiki.chasmsoftware.org/index.php/SoftwareLicense. Commercial users should contact the Johns Hopkins Technology Transfer office. | Yes | Yes | Yes | DbNSFP VEP Plugin
phastCons100way_vertebrate | phastCons conservation score based on the multiple alignments of 100 vertebrate genomes (including human). The larger the score, the more conserved the site. | Yes | Yes | Yes | DbNSFP VEP Plugin
phastCons100way_vertebrate_rankscore | phastCons100way_vertebrate scores were ranked among all phastCons100way_vertebrate scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of phastCons100way_vertebrate scores in dbNSFP. | Yes | Yes | Yes | DbNSFP VEP Plugin
phastCons46way_placental | phastCons conservation score based on the multiple alignments of 33 placental mammal genomes (including human). The larger the score, the more conserved the site. | Yes | Yes | Yes | DbNSFP VEP Plugin
phastCons46way_placental_rankscore | phastCons46way_placental scores were ranked among all phastCons46way_placental scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of phastCons46way_placental scores in dbNSFP. | Yes | Yes | Yes | DbNSFP VEP Plugin
phastCons46way_primate | phastCons conservation score based on the multiple alignments of 10 primate genomes (including human). The larger the score, the more conserved the site. | Yes | Yes | Yes | DbNSFP VEP Plugin
phastCons46way_primate_rankscore | phastCons46way_primate scores were ranked among all phastCons46way_primate scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of phastCons46way_primate scores in dbNSFP. | Yes | Yes | Yes | DbNSFP VEP Plugin
phyloP100way_vertebrate | phyloP (phylogenetic p-values) conservation score based on the multiple alignments of 100 vertebrate genomes (including human). The larger the score, the more conserved the site. | Yes | Yes | Yes | DbNSFP VEP Plugin
phyloP100way_vertebrate_rankscore | phyloP100way_vertebrate scores were ranked among all phyloP100way_vertebrate scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of phyloP100way_vertebrate scores in dbNSFP. | Yes | Yes | Yes | DbNSFP VEP Plugin
phyloP46way_placental | phyloP (phylogenetic p-values) conservation score based on the multiple alignments of 33 placental mammal genomes (including human). The larger the score, the more conserved the site. | Yes | Yes | Yes | DbNSFP VEP Plugin
phyloP46way_placental_rankscore | phyloP46way_placental scores were ranked among all phyloP46way_placental scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of phyloP46way_placental scores in dbNSFP. | Yes | Yes | Yes | DbNSFP VEP Plugin
phyloP46way_primate | phyloP (phylogenetic p-values) conservation score based on the multiple alignments of 10 primate genomes (including human). The larger the score, the more conserved the site. | Yes | Yes | Yes | DbNSFP VEP Plugin
phyloP46way_primate_rankscore | phyloP46way_primate scores were ranked among all phyloP46way_primate scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of phyloP46way_primate scores in dbNSFP. | Yes | Yes | Yes | DbNSFP VEP Plugin
Condel | Consensus deleteriousness score for an amino acid substitution based on SIFT and PolyPhen-2 | Yes | Yes | Yes | Condel VEP Plugin
Grantham | Grantham Matrix score - Grantham, R. Amino Acid Difference Formula to Help Explain Protein Evolution, Science 1974 Sep 6;185(4154):862-4. | Yes | Yes | Yes | Grantham VEP Plugin

## Provenance Report
The provenance report is designed to ensure reproducibility of processing by documenting the details of the input and
output files and the versions of software tools applied to them. The provenance report is named 
`<pipeline run ID>_<sample name>.provenance.pdf`, and documents the following:

* The time and date of the beginning and end of processing the sample
* The patient identifier
* The timestamps and file sizes of
    * The input FASTQ files
    * The initial alignment file (BAM file)
    * The final alignment file (ending in .recal.bam)
    * The raw VCF file
    * The annotation file produced by Annovar (CSV format)
    * The final Excel report containing variant calls
* The version numbers of all the critical software components used in generating the results

## Gap Report

The gap report shows additional annotated information regarding overlapping coding regions for each sample. This report
 is entitled `<pipeline run ID>_<sample name>.gap.csv` in the `results` directory.
 
The gap report has the following fields:

Field Name | Description | Source
--- | --- | ---
Chr | chromosome of gap | bed file
Gene | gene of gap | bed file
Start | start position of gap | calculated based on first position of gap with coverage below threshold (inclusive)
End | end position of gap | last position of gap with coverage below threshold (exclusive)
Min Cov | Lowest coverage in gap | calculated across gap
Max Cov | Highest coverage in gap | calculated across gap (must be less than or equal to max_low_coverage)
Median Cov | Middle value of coverage | calculated across gap
Mean Cov | Average coverage across gap | calculated across gap
Width | Width of gap in bases | calculated across gap
Tx Name | Transcript name in overlapping coding region | refgene
Strand | Strand of overlapping region | refgene
CDS Distance | Distance to nearest CDS region (-ve if direction is not same as strand) | refgene
CDS Overlap Start | Start of overlap relative to exon (CDS) starting from strang | refgene
CDS Overlap End | End of overlap (inclusive) relative to exon (CDS) | refgene
CDS Segment Start | Start of overlap relative to CDS start, only including exons | refgene
CDS Segment End | End of overlap relative to CDS start, only including exons | refgene
AA Overlap Start | Codon position start of overlap | refgene
AA Overlap End | Codon position end of overlap | refgene
Exon Number | Exon number | refgene
Exon Rank | Exon rank in the transcript counting in the direction of the strand | Same as exon number for + strand, else counting starts from end.
(Bed filename) | List of overlap names separated by ; |  

