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

// all the core pipeline stages in the pipeline
load 'pipeline_stage_initialize.groovy'
load 'pipeline_stage_alignment.groovy'
load 'pipeline_stage_variant_calling.groovy'
load 'pipeline_stage_annotation.groovy'
load 'pipeline_stage_reports.groovy'

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
    initialize_batch_run +

    // for each analysis profile we run the main pipeline in parallel
    ANALYSIS_PROFILES * 
    [
        initialize_profiles +
        
        samples * // for each sample...
        [
            // phase 1. data pre-processing for each sample: alignment, mark duplicates, indel realignment, base recalibration -> analysis ready reads
            align_sample + 
            [ 
                // phase 2. variant calling
                variant_discovery + 

                // phase 3. annotation
                variant_annotation + 

                // phase 4. sample specific reports
                sample_reports, sample_reports_extra
            ] +
            sample_checks
        ] + 
        qc_excel_report
   ] +

   // produce the output spreadsheet, 1 per analysis profile
   ANALYSIS_PROFILES * 
   [ 
       set_target_info +  
       [ 
           vcf_to_excel, 
           family_vcf 
       ]
   ] +

   // Produce a mini bam for each variant to help investigate individual variants
   samples * 
   [ 
       variant_bams, 
       filtered_on_exons + 
       index_bam 
   ] +

   // And then finally write the provenance report (1 per sample)
   samples * 
   [ 
       provenance_report
   ] +
   
   // clean up
   finish_batch_run
}

