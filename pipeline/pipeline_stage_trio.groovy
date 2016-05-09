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
    def additional_samples = sample_info[sample].pedigree.tokenize(';')[0].tokenize('=')[1].tokenize(',');
    // from_list = add
    def from_list = additional_samples.collect { "variants/${it}.combined.g.vcf" }
    from_list.add("variants/${sample}.combined.g.vcf")
    stage_status("trio_analysis_phase_2", "sample: ${sample}; trio samples: ${additional_samples}; from list: ${from_list}", sample)

    output.dir="variants"

    produce ("${sample}.trio.genotype.raw.vcf") {
        from(from_list) {
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
        }
    } // produce
    stage_status("trio_analysis_phase_2", "exit", sample)
}

genotype_refinement_trio = {
    doc "following the gatk workflow, these steps improve the accuracy of trio based calls"
    stage_status("genotype_refinement_trio", "enter", sample)
    // more details at https://www.broadinstitute.org/gatk/guide/article?id=4723, https://www.broadinstitute.org/gatk/guide/article?id=4727

    // Step 1: Derive posterior probabilities of genotypes
    // java -Xmx20g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T CalculateGenotypePosteriors -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --supporting ref_files/1000G_phase3_v4_20130502.sites.hg19.vcf -ped txxxx.ped -V txxxx.genotype.raw.vcf -o txxxx.genotype.raw.postCGP.vcf
    // Step 2: Filter low quality genotypes
    // java -Xmx10g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T VariantFiltration -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta -V txxxx.genotype.raw.postCGP.vcf -G_filter "GQ < 20.0" -G_filterName lowGQ -o txxxx.genotype.raw.postCGP.GQfilter.vcf
    // Step 3: Annotate possible de novo mutations
    // java -Xmx10g -jar /usr/local/gatk/4.5/GenomeAnalysisTK.jar -T VariantAnnotator -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta -V txxxx.genotype.raw.postCGP.GQfilter.vcf -A PossibleDeNovo -ped txxxx.ped -o txxxx.genotype.raw.postCGP.GQfilter.deNovos.vcf
    output.dir="variants"
    def safe_tmp_dir = [TMPDIR, UUID.randomUUID().toString()].join( File.separator )
    produce("${sample}.${analysis}.refined.vcf") {
        stage_status("genotype_refinement_trio", "sources are variants/${sample}.${analysis}.genotype.raw.vcf and results/${run_id}_family_${sample}.ped", sample)
        from("variants/${sample}.${analysis}.genotype.raw.vcf", "results/${run_id}_family_${sample}.ped") {
            exec """
                mkdir -p "$safe_tmp_dir"

                java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T CalculateGenotypePosteriors -R $REF --supporting $TRIO_REFINEMENT_SUPPORTING -ped ${input.ped} -V ${input.vcf} -o "${safe_tmp_dir}/postCGP.vcf"

                java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $REF -V "${safe_tmp_dir}/postCGP.vcf" -G_filter "GQ < 20.0" -G_filterName lowGQ -o "${safe_tmp_dir}/GQfilter.vcf"

                java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -T VariantAnnotator -R $REF -V "${safe_tmp_dir}/GQfilter.vcf" -A PossibleDeNovo -ped ${input.ped} -o ${output}

                rm -r "$safe_tmp_dir"
            """, "genotype_refinement"
        }
    }

    stage_status("genotype_refinement_trio", "exit", sample)
}

