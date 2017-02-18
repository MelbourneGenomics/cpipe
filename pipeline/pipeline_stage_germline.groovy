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
load 'pipeline_helpers.groovy'

/////////////////////////////////////////////////////////////////////////////////
germline_analysis_phase_1 = segment {
    // all samples do this
    variant_discovery // stage_variant_calling, result is sample.combined.g.vcf
}

/////////////////////////////////////////////////////////////////////////////////
// given the g.vcf of the individual, convert this to a genotype.raw.vcf
// we do this for probands and individuals, but not samples that are just for trios
germline_genotype_gvcfs = {
    stage_status('germline_genotype_gvcfs', 'enter', sample);
    // joint_call_individual
    // java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --disable_auto_index_creation_and_locking_when_reading_rods --num_threads 1 --variant 00NA12877.hap.raw.g.vcf --variant 00NA12878.hap.raw.g.vcf --variant 00NA12879.hap.raw.g.vcf --out txxxx.genotype.raw.vcf -ped txxxx.ped -log txxxx.GenotypeGVCFs.log --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf -G Standard -A AlleleBalance -A AlleleBalanceBySample -A DepthPerAlleleBySample -A GCContent -A GenotypeSummaries -A HardyWeinberg -A LikelihoodRankSumTest -A MappingQualityZero -A SampleList -A SpanningDeletions -A StrandBiasBySample -A TandemRepeatAnnotator -A VariantType -A TransmissionDisequilibriumTest
    // 23-may-2016
    // java -Xmx24g -jar /usr/local/gatk/3.5/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /vlsci/VR0320/shared/production/1.0.4/hg19/ucsc.hg19.fasta --disable_auto_index_creation_and_locking_when_reading_rods --num_threads 1 --variant 00NA12877.hap.raw.g.vcf --variant 00NA12878.hap.raw.g.vcf --variant 00NA12879.hap.raw.g.vcf --out txxxx.genotype.raw.fewerA.vcf -ped txxxx.ped -log txxxx.GenotypeGVCFs.fewerA.log --dbsnp /vlsci/VR0320/shared/production/1.0.4/hg19/dbsnp_138.hg19.vcf -A AlleleBalance -A VariantType
    output.dir="variants"

    var call_conf:5.0

    produce("${sample}.individual.genotype.raw.vcf") {
        from("variants/${sample}.hc.g.vcf") {
            exec """
                $JAVA -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs
                    -R $REF
                    --disable_auto_index_creation_and_locking_when_reading_rods
                    --num_threads $threads
                    --variant $input
                    --out $output
                    --logging_level INFO
                    --dbsnp $DBSNP
                    -A AlleleBalance
                    -A VariantType
                    -stand_call_conf $call_conf 
            """, "gatk_genotype"
        }
    }
    stage_status('germline_genotype_gvcfs', 'exit', sample);
}

/////////////////////////////////////////////////////////////////////////////////
// this is the variant only pathway
// we already have the vcf, just forward it on to the next stage
germline_genotype_vcfs = {
    stage_status('germline_genotype_vcfs', 'enter', sample);
    output.dir="variants"
    produce("${sample}.individual.genotype.raw.vcf") {
        from("variants/${sample}.hc.vcf") {
           exec """
               cp "$input" "$output"
           """
        }
    }

    stage_status('germline_genotype_vcfs', 'exit', sample);
}

/////////////////////////////////////////////////////////////////////////////////
genotype_refinement_individual = {
    doc "skips the refinement steps and passes on the vcf to the next stage unchanged"
    output.dir="variants"
    produce("${sample}.${analysis}.refined.vcf") {
        from("variants/${sample}.${analysis}.combined.genotype.vcf") {
            exec """
                cp $input $output
            """
        }
    }
}

if (GATK_VARIANT_ONLY) {
    germline_analysis_phase_2 = segment {
        germline_genotype_vcfs +
        filter_variants +
        merge_variants
    }
}
else {
    germline_analysis_phase_2 = segment {
        germline_genotype_gvcfs +
        filter_variants +
        merge_variants
    }
}
