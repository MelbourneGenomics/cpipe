// vim: ts=4:sw=4:expandtab:cindent
/////////////////////////////////////////////////////////////////////////////////
//
// This file is part of Cpipe.
// 
// Cpipe is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, under version 3 of the License, subject
// to additional terms compatible with the GNU General Public License version 3,
// specified in the LICENSE file that is part of the Cpipe distribution.
//
// Cpipe is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Cpipe.  If not, see <http://www.gnu.org/licenses/>.
//
/////////////////////////////////////////////////////////////////////////////////
//
// Cpipe Main Pipeline Script
//
/////////////////////////////////////////////////////////////////////////////////

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
// load 'haloplex.groovy'

sample_metadata_file = correct_sample_metadata_file( args[0] ) // fix syntax issues and update sample_metadata_file

try {
  sample_info = SampleInfo.parse_mg_sample_info(sample_metadata_file)
}
catch (RuntimeException e) {
  sample_info = SampleInfo.parse_sample_info(sample_metadata_file)
}

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

run {
    // Check the basic sample information first
    check_sample_info +  // check that fastq files are present
    check_tools +
    update_gene_lists + // build new gene lists by adding sample specific genes to cohort

    // Create a single BED that contains all the regions we want to call
    // variants in
    create_combined_target + 
    create_synonymous_target + // regions where synonymous snvs are not filtered
    build_capture_stats + // how well covered genes are by the capture

    generate_pipeline_id + // make a new pipeline run ID file if required

    // For each analysis profile we run the main pipeline in parallel
    ANALYSIS_PROFILES * [

        set_target_info + 

        init_analysis_profile +

        create_splice_site_bed +
        
        // The first phase is to perform alignment and variant calling for each sample
        samples * [
               set_sample_info +
                   "%.gz" * [ fastqc ] + check_fastqc +
                   ~"(.*)_R[0-9][_.].*fastq.gz" * [ trim_fastq + align_bwa + index_bam + cleanup_trim_fastq ] +
                   merge_bams +
                   dedup + 
                   cleanup_initial_bams +
                   realignIntervals + realign + index_bam +
                   bsqr_recalibration + index_bam +
				   cleanup_intermediate_bams +
                       [
                         call_variants_gatk + call_pgx + merge_pgx +
                         filter_variants + merge_variants +
                         annotate_vep + index_vcf +
                         annovar_table +
                         [ 
                             add_to_database, 
                             augment_condel + annotate_significance
                         ]  +
                         [ calc_coverage_stats + check_ontarget_perc, calculate_fragment_statistics ] + [ summary_report, exon_qc_report, gap_report ],
                         gatk_depth_of_coverage,
                         insert_size_metrics
                       ]
                   + check_coverage
                   + check_karyotype
        ] + qc_excel_report 
   ] + 

   // The 3rd phase is to produce the output spreadsheet, 1 per analysis profile
   ANALYSIS_PROFILES * [ set_target_info +  [ vcf_to_excel, family_vcf ] ] +

   // Produce a mini bam for each variant to help investigate individual variants
   samples * [ variant_bams ] +

   // And then finally write the provenance report (1 per sample)
   samples * [ provenance_report /* , annovar_to_lovd */ ] +
   
   // And report on similarity between samples
   sample_similarity_report +

   // check overall quality of results
   validate_batch +

   // update metadata and pipeline ID
   create_sample_metadata
}

