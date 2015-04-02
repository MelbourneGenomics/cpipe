// vim: ts=4:sw=4:expandtab:cindent:number
////////////////////////////////////////////////////////////
//
// Melbourne Genomics Variant Calling Pipeline
//
// Pipeline stage definitions for exome variant calling 
// pipeline. See pipeline.groovy for more information.
//
// Author: Simon Sadedin, MCRI
//         Members of Melbourne Genomics
//
// Copyright Melbourne Genomics Health Alliance members. All rights reserved.
//
// DISTRIBUTION:
//
// This source code should not be distributed to a third party without prior
// approval of the Melbourne Genomics Health Alliance steering committee (via
// Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
/// 
////////////////////////////////////////////////////////

add_splice_variants = {
    doc "Annovar seems to leave out lot of variants that VCGS would consider interesting splice variants, so we add them in afterwards."
    output.dir = "variants"

    if(!splice_region_window)
        return

    msg "Adding splice variants from genome summary to exome summary ($input.csv) ..."

    filter("addsplice") {
        from("exome_summary.csv","genome_summary.csv") {
            exec "cp $input1 $output; python $SCRIPTS/add_splice_variants.py $exon_bed_file $input2 $splice_region_window >> $output"
        }
    }
}


@filter("aug")
augment_transcript_ids = {
    doc "Add an additional column indicating the VCGS transcript identifier for the VCGS transcript / isoform affected by each variant"
    output.dir="variants"

    msg "Augmenting Annovar output with extra columns ($input.csv) ..."
    exec "python $SCRIPTS/augment_transcripts.py $transcripts_file $input.csv $input.exonic_variant_function > $output.csv"
}


multiple_annovar = segment {
    [
        annovar_summarize.using(source:"refgene") + add_splice_variants + augment_transcript_ids,
        annovar_summarize.using(source:"knowngene") + add_splice_variants + augment_transcript_ids
    ] + merge_annovar_reports 
}



