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


//////////////////////////////////////////////////////////////////////
// stages
//////////////////////////////////////////////////////////////////////

// at this point all samples have .g.vcf
// want to generate the trio based genotype.raw.vcf
// assumes this is called for all proband samples
trio_analysis_phase_2 = {
    stage_status("trio_analysis_phase_2", "enter", sample)
    // extract parent sample ids
    additional_samples = sample_info[sample].pedigree.tokenize(';')[0].tokenize('=')[1].tokenize(',');
    // from_list = add
    from_list = additional_samples.collect { "${it}.combined.g.vcf" }
    from_list.add("${sample}.combined.g.vcf")
    stage_status("trio_analysis_phase_2", "trio samples: ${additional_samples}; from list: ${from_list}", sample)
    println("from list: $from_list")

    output.dir="variants"

    from(from_list) produce ("${sample}.trio.genotype.raw.vcf") {
        stage_status("trio_analysis_phase_2", "inputs: ${inputs}", sample)
        additional_variant_params = inputs.collect { "--variant $it" }.join(' ')
        exec """
            java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs
                -R $REF
                --disable_auto_index_creation_and_locking_when_reading_rods
                --num_threads $threads
                $additional_variant_params
                --out $output
                --logging_level INFO
                --dbsnp $DBSNP
                -G Standard
                -A AlleleBalance
                -A AlleleBalanceBySample
                -A DepthPerAlleleBySample
                -A GCContent
                -A GenotypeSummaries
                -A HardyWeinberg
                -A LikelihoodRankSumTest
                -A MappingQualityZero
                -A SampleList
                -A SpanningDeletions
                -A StrandBiasBySample
                -A TandemRepeatAnnotator
                -A VariantType
                -A TransmissionDisequilibriumTest
        """, "gatk_genotype"
    } // produce
    stage_status("trio_analysis_phase_2", "exit", sample)
}
