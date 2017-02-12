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

// remove spaces from gene lists and point to a new sample metadata file
// note that this isn't run through bpipe
correct_sample_metadata_file = { sample_meta_data_file ->
    def target = new File('results')
    if( !target.exists() ) {
        target.mkdirs()
    }
    
    def outputFile = "results/samples.corrected"
    
    if(file(outputFile).exists() && file(outputFile).lastModified() > file(sample_meta_data_file).lastModified()) {
        println "Skip creating corrected meta data file (up-to-date)"
    }
    else {
        println "Creating corrected meta data file"
        println([ "sh", "-c", "python $SCRIPTS/correct_sample_metadata_file.py < $sample_meta_data_file > $outputFile" ].execute().text)
    }
    return outputFile
}

/////////////////////////////////////////////////////////
// helper functions
/////////////////////////////////////////////////////////

// classify samples as singleton, trio, or member of trio
List find_sample_types(sample_info) {
    println "finding samples"

    // all samples
    all_samples = sample_info.keySet()

    // proband samples
    trio_samples = []
    proband_samples = all_samples.findAll { 
        if (sample_info[it].pedigree != "" && sample_info[it].pedigree != "exclude" && sample_info[it].pedigree != "import") {
            new_members = sample_info[it].pedigree.tokenize(';')[0].tokenize('=')[1].tokenize(','); // fid=na12877,na12878
            trio_samples.addAll(new_members); // add to members
            return true;
        }
        else {
            return false;
        }
    }

    // only analyze non-trio
    // individual_samples = all_samples.collect()
    // individual_samples.removeAll(trio_samples) // includes probands

    // analyze unless exclude
    individual_samples = all_samples.findAll { 
        if (sample_info[it].pedigree == "exclude") {
            return false; // exclude
        }
        else {
            return true; // include
        }
    }
    println "done finding samples"
    return [ all_samples, proband_samples, trio_samples, individual_samples ]
}

// log changes to stage status
void stage_status(stage_name, stage_status, sample) {
    String current = new Date().format("yyMMdd-HHmmss")
    println("${current}: ${stage_name}: ${stage_status} (${sample})")
}

/////////////////////////////////////////////////////////
// common stages
/////////////////////////////////////////////////////////

filter_variants = {
    doc "Select only variants in the genomic regions defined for the $target_name target"
    output.dir="variants"

    stage_status("filter_variants", "enter", sample)

    def pgx_flag = ""
    if(file("../design/${target_name}.pgx.vcf").exists()) {
        pgx_flag = "-L ../design/${target_name}.pgx.vcf"
    }

    msg "Filtering variants - finding INDELs"
    exec """
        java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar 
             -R $REF
             -T SelectVariants 
             --variant $input.vcf 
             -L $target_bed_file.${sample}.bed $pgx_flag --interval_padding $INTERVAL_PADDING_SNV
             --selectTypeToInclude SNP --selectTypeToInclude MIXED --selectTypeToInclude MNP --selectTypeToInclude SYMBOLIC --selectTypeToInclude NO_VARIATION
             -o $output.snv
    """

    msg "Filtering variants - finding SNVs"
    exec """
        java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar 
             -R $REF
             -T SelectVariants 
             --variant $input.vcf 
             -L $target_bed_file.${sample}.bed $pgx_flag --interval_padding $INTERVAL_PADDING_INDEL
             --selectTypeToInclude INDEL
             -o $output.indel
    """
    stage_status("filter_variants", "exit", sample)
}

merge_variants_gvcf = {
    doc "Merge SNVs and INDELs"
    output.dir="variants"
    stage_status("merge_variants_gvcf", "enter", sample)

    produce("${sample}.combined.g.vcf") {
        exec """
            java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T CombineGVCFs
            -R $REF
            --variant:indel $input.indel
            --variant:snv $input.snv
            --out $output.combined.g.vcf
        """
    }
    stage_status("merge_variants_gvcf", "exit", sample)
}


merge_variants = {
    doc "Merge SNVs and INDELs"
    output.dir="variants"
    stage_status("merge_variants", "enter", "${sample} ${branch.analysis}")

    msg "Merging SNVs and INDELs"
    produce("${sample}.${analysis}.combined.genotype.vcf") {
        exec """
            java -Xmx3g -jar $GATK/GenomeAnalysisTK.jar
            -T CombineVariants
            -R $REF
            --variant:indel $input.indel
            --variant:snv $input.snv
            --out $output.vcf
            --setKey set
            --genotypemergeoption UNSORTED
         """
    }
    stage_status("merge_variants", "exit", "${sample} ${branch.analysis}")
}

