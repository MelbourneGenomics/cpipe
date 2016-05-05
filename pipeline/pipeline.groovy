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
load 'pipeline_helpers.groovy'

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

/////////////////////////////////////////////////////////
sample_metadata_file = correct_sample_metadata_file(args[0]) // fix syntax issues and update sample_metadata_file

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


// trio members
(all_samples, proband_samples, trio_samples, individual_samples) = find_sample_types(sample_info)
println "all samples: ${all_samples}; proband_samples: ${proband_samples}; trio_samples: ${trio_samples}; individual_samples: ${individual_samples}"

// due to bpipe not allowing empty branch lists https://github.com/ssadedin/bpipe/issues/180
// we have to explicitly mark a list as empty
// proband is the only one that could feasibly be empty and is used in a branch list
EMPTY_MARKER = "** this is not a real sample **"
if (proband_samples.size() == 0) {
    proband_samples.add(EMPTY_MARKER)
}

// all the core pipeline stages in the pipeline
load 'pipeline_stage_initialize.groovy' // preparation of a batch
load 'pipeline_stage_alignment.groovy' // generate a bam
load 'pipeline_stage_variant_calling.groovy' // find variants
load 'pipeline_stage_variant_analysis.groovy' // filter, normalize, annotate, post process
load 'pipeline_stage_reports.groovy'

// specific to type of analysis
load 'pipeline_stage_germline.groovy'
load 'pipeline_stage_trio.groovy'

set_sample_name_without_target = {
    if (branch.name == EMPTY_MARKER) {
        stage_status('set_sample_name', 'skipping empty branch', branch.name)
        succeed "This is a dummy branch. Not an error."
    }
    branch.sample = branch.name
}

set_sample_name = {
    if (branch.name == EMPTY_MARKER) {
        stage_status('set_sample_name', 'skipping empty branch', branch.name)
        succeed "This is a dummy branch. Not an error."
    }
    branch.sample = branch.name

    // terminate the branch if the profile doesn't match
    if(sample_info[branch.sample].target != branch.target_name) {
        // This is expected because every file is processed for every target/flagship
        succeed "skipping sample $sample for target $target_name"
    }
}

set_analysis_type = {
    println("updating analysis from ${branch.analysis} to ${new_analysis}")
    branch.analysis = new_analysis
}

set_analysis_type_individual = {
    println("updating analysis from ${branch.analysis} to individual")
    branch.analysis = "individual"
}

set_analysis_type_trio = {
    println("updating analysis from ${branch.analysis} to trio")
    branch.analysis = "trio"
}

run {

    initialize_batch_run + // some overall checks, overall target region, ped files, pipeline run ID

    // for each analysis profile we run the main pipeline in parallel
    ANALYSIS_PROFILES * 
    [
        initialize_profiles + // setup target regions
        
        all_samples * // for each sample...
        [
            // --- module 1. data pre-processing for each sample: ---
            // the goal of this module is an analysis ready BAM
            // alignment, mark duplicates, indel realignment, base recalibration -> analysis ready reads
            // loosely based on gatk workflow
            set_sample_name +
            analysis_ready_reads + // pipeline_stages_alignment

            // --- module 2. variant discovery
            // the goal of this module is a raw VCF
            // each sample is passed onto each type of enabled analysis and processed if relevant
            germline_analysis_phase_1 + // haplotypecaller without a ped for all samples (variant_discovery). Result -> all samples have .g.vcf

            // generate reports and do checks based on bam alignment
            analysis_ready_reports +
            analysis_ready_checks // reports
        ] +

        // each type of analysis can do a second phase after all phase 1 stages have finished
        // i.e. because trio analysis is dependent on germline analysis phase 1
        // result is individual.genotype.raw.vcf and trio.genotype.raw.vcf
        [
            proband_samples *
            [
                set_sample_name + trio_analysis_phase_2 // .trio.genotype.raw.vcf
            ],
            individual_samples * // individuals and probands
            [
                set_sample_name + germline_analysis_phase_2 // .individual.genotype.raw.vcf
            ]
        ] +

        // --- module 3. fitering and annotation
        [ 
            proband_samples *
            [
                // set_sample_name + set_analysis_type.using(new_analysis: "trio") + variant_analysis
                set_sample_name + set_analysis_type_trio + genotype_refinement_trio + variant_analysis
            ],
            individual_samples * // individuals and probands
            [
                // set_sample_name + set_analysis_type.using(new_analysis: "individual") + variant_analysis
                set_sample_name + set_analysis_type_individual + genotype_refinement_individual + variant_analysis
            ]
        ]

        // qc_excel_report deprecated, see analysis_ready_reports
   ] +

   // --- module 4. post-processing
   // Produce a mini bam for each variant to help investigate individual variants
   // * mini-bams are no longer produced

   // And then finally write the provenance report (1 per sample)
   all_samples *
   [ 
       provenance_report
   ] +

   // clean up
   finish_batch_run
}

