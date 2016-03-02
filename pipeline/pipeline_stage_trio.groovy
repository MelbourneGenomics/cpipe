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

call_variants_child_trio = {
    doc "Call SNPs/SNVs using GATK Haplotype Caller for all samples"
    output.dir="variants"

    // Default values of confidence thresholds
    // come from the Broad web site. However
    // these may be higher than suitable in our context
    var call_conf:5.0, 
        emit_conf:5.0

    transform("bam") to("vcf") {
        exec """

            $JAVA -Xmx6g -jar $GATK/GenomeAnalysisTK.jar -T HaplotypeCaller
                   -R $REF 
                   -I $input.bam 
                   -L $COMBINED_TARGET $splice_region_bed_flag
                   --interval_padding $INTERVAL_PADDING_CALL
                   --dbsnp $DBSNP 
                   --emitRefConfidence GVCF
                   -stand_call_conf $call_conf 
                   -stand_emit_conf $emit_conf
                   -dcov 1600 
                   -l INFO 
                   -A AlleleBalance 
                   -A Coverage 
                   -A FisherStrand 
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}


//////////////////////////////////////////////////////////////////////
// segments
//////////////////////////////////////////////////////////////////////

trio_analysis_phase_1 = segment {
    dummy
}

trio_analysis_phase_2 = segment {
    // TODO if (sample.is_trio_child()) {
    //    call_variants_child_trio +
    //    joint_call + 
    //    variant_annotation 
    // TODO }
    dummy
}
