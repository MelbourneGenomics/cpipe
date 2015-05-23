Cpipe Designs (Analysis Profiles)
----------------------------------

Cpipe can conduct many different kinds of analysis. The configuration for each
analysis is called an "Analysis Profile", and is stored in a dedicated subdirectory
of this folder. Adding new analysis profiles should be performed by running
the 'pipeline/scripts/new_target_region.sh' script, which performs a number of checks
and some pre-annotation steps that are required before Cpipe runs on a data set.

# Analysis Profile Structure

## Diagnostic Target Region
The key component of an analysis profile is the set of genes (or genomic regions) to 
be analysed. These are the genes for which variant calls will be written into
the final output file. These regions (or genes) are called the "diagnostic
target region", and are separate to the "exome target region". To conduct an
analysis, both the diagnostic region and the exome regions must be defined.

The diagnostic target region is supplied as a BED file that defines the regions
you wish to analyse.

## Prioritised Genes
An analysis profile can define tiered gene categories for the genes within the
diagnostic target region. These gene categories are defined in a file called

  PROFILE_ID.genes.txt

where each line has an HGNC gene symbol and an integer priority, separated by a
tab character. Higher numbers are treated as more important. By convention, the
following interpretation is used:

    0 - Exclude from reporting (not yet implemented!)
    1 - Default
    2 - Genes with known clinical diagnostic variants (across different diseases)

Category 3 and 4 are by convention specified per-patient in the samples.txt 
file and would not be entered into the analysis profile gene category file.

  3 - Genes implicated as likely to be significant for the patient (eg: by 
      phenotype prediction software or implicated by clinical diagnosis for
      patient)
  4 - Genes specifically identified by a clinician based on experience or 
      clinical acumen

## Analysis Specific Settings

Cpipe offers many different adjustable settings and thresholds. For example, the 
minimum coverage depth for processing samples, the maximum PCR duplication rate,
the width of window out from exons used to identify splice variants, etc. In
addition, entire pipeline stages and tools can be replaced on a per analysis
profile basis. All of these settings and customisations are performed by the
"settings file" which is named according to the following convention:

  PROFILE_ID.settings.txt


# Creating a Batch

Creating an analysis batch is performed using the
"pipeline/scripts/create_batch.sh" script. This script requires specification of
the analysis profile to use for the analysis, as well as the exome (or panel)
target region to use for the analysis. The analysis profile should be specified
by its profile ID. The exome target region can be specified in one of the
following ways:

  * an absolute path to a BED file defining the target regions
  * a relative path to a BED file inside the Cpipe directory
  * an unqualified path to a BED file that is inside the analysis profile
    directory

Further instructions on creating the batches are printed out by the
create_batch.sh script itself.
