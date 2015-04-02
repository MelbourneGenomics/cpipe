// vim: ts=4:sw=4:expandtab:cindent
////////////////////////////////////////////////////////////////////////////
//
// Melbourne Genomics Variant Calling Pipeline
//
// This pipeline executes the standard variant calling analysis for 
// analysing data for the MGHA project.
//
// The best available documentation for setting up the pipeline is located
// here:
//
//   https://sites.google.com/site/melbournegenomics/how-tos/setting-up-the-pipeline
// 
// Usage:
//
//   bpipe run ../../../pipeline/pipeline.groovy ../samples.txt
// 
// Author: Simon Sadedin, MCRI
//         Members of the Melbourne Genomics
// 
// Copyright Melbourne Genomics Health Alliance members. All rights reserved.
//
// DISTRIBUTION:
//
// This source code should not be distributed to a third party without prior
// approval of the Melbourne Genomics Health Alliance steering committee (via
// Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
//
////////////////////////////////////////////////////////////////////////////

about title: "Melbourne Genomics Demonstration Project Pipeline"

// Load the default configuration
load 'config.groovy'

// Local file can set EXOME_TARGET and ANALYSIS_PROFILES
if(file("../target_regions.txt").exists())  {
    load '../target_regions.txt'
}

requires EXOME_TARGET : """
        The exome target regions. This should be the whole regions targeted
        for capture by the exome kit. Note that the regions for analysis
        may be a subset, but you should always specify the whole exome
        region here.
    """

// All the core pipeline stages in the pipeline
load 'pipeline_stages_config.groovy'

// VCGS specific stages
load 'vcgs.groovy'

sample_metadata_file = args[0]
sample_info = SampleInfo.parse_sample_info(args[0])

// We are specifying that each analysis takes place inside a fixed file structure
// where the parent directory is named according to the batch name. Thus we
// can infer the batch name from the name of the parent directory.
// 
// Note: this variable can be overridden by passing a parameter to bpipe in case
// you are running in a different location.
batch = new File("..").canonicalFile.name

// Extract the analysis profiles from the sample information
ANALYSIS_PROFILES = sample_info*.value*.target as Set

samples = sample_info.keySet()

annovar = segment {
    annovar_summarize + add_splice_variants
}

run {
    // Check the basic sample information first
    check_sample_info + check_tools +

    // Create a single BED that contains all the regions we want to call
    // variants in
    create_combined_target + 

    // For each analysis profile we run the main pipeline in parallel
    ANALYSIS_PROFILES * [

        set_target_info +
        
        // The first phase is to perform alignment and variant calling for each sample
        samples * [
               set_sample_info +
                   "%.gz" * [ fastqc ] + check_fastqc +
                   ~"(.*)_R[0-9][_.].*fastq.gz" * [ trim_fastq + align_bwa + index_bam + cleanup_trim_fastq ] +
                   merge_bams +
                   dedup + index_bam + 
                   cleanup_initial_bams +
                   realignIntervals + realign + index_bam +
                   recal_count + recal + index_bam +
				   cleanup_intermediate_bams +
                       [
                         call_variants_hc + call_pgx + merge_pgx +
                            filter_variants + 
                            annotate_vep + index_vcf +
                            annovar +
                            [ 
                               add_to_database, 
                               augment_condel + annotate_significance, 
                               calculate_cadd_scores
                            ] + augment_cadd +
                         calc_coverage_stats + summary_pdf, 
                         gatk_depth_of_coverage,
                         insert_size_metrics
                       ]
                   + check_coverage
                   + check_karyotype
        ] + qc_excel_report
   ] + 

   // The 3rd phase is to produce the output spreadsheet, 1 per analysis profile
   ANALYSIS_PROFILES * [ set_target_info +  [ vcf_to_excel, family_vcf ] ] +

   // And then finally write the provenance report (1 per sample)
   samples * [ provenance_report /* , annovar_to_lovd */ ] +
   
   // And report on similarity between samples
   sample_similarity_report
}

