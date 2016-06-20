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
load "pipeline_helpers.groovy"

///////////////////////////////////////////////////////////////////
// stages
///////////////////////////////////////////////////////////////////
calc_coverage_stats = {
    doc "Calculate coverage across a target region using Bedtools"
    stage_status("calc_coverage_stats", "enter", sample)
    output.dir="qc"

    var MIN_ONTARGET_PERCENTAGE : 50

    // transform("bam") to([sample + ".cov.gz", sample + ".exome.gz", sample + ".ontarget.txt"]) {
    produce("${sample}.cov.gz", "${sample}.exome.gz", "${sample}.ontarget.txt") {
        // only calculate coverage for bases overlapping the capture
        def safe_tmp_dir = [TMPDIR, UUID.randomUUID().toString()].join( File.separator )

        // coverage calculation for qc_report.py
        // 1 determine the intersection between the target region and the exome
        // 2 calculate the coverage in the bam file over this intersection
        // 3 calculate the coverage in the bam file over the exome
        // 4 calculate reads on target over the combined target region
        exec """
          mkdir -p "$safe_tmp_dir"
        
          $BEDTOOLS/bin/bedtools intersect -a $target_bed_file.${sample}.bed -b $EXOME_TARGET | sort -k 1,1 -k2,2n > "$safe_tmp_dir/intersect.bed"

          $BEDTOOLS/bin/bedtools bamtobed -i $input.recal.bam | sort -k 1,1 -k2,2n > "$safe_tmp_dir/bam.bed"

          $BEDTOOLS/bin/coverageBed -d -sorted -a "$safe_tmp_dir/intersect.bed" -b "$safe_tmp_dir/bam.bed" | gzip > $output.cov.gz

          sort -k 1,1 -k2,2n < $EXOME_TARGET > "$safe_tmp_dir/exome.bed"

          $BEDTOOLS/bin/coverageBed -d -sorted -a "$safe_tmp_dir/exome.bed" -b "$safe_tmp_dir/bam.bed" | gzip > $output2.exome.gz
        
          $SAMTOOLS/samtools view -L $COMBINED_TARGET $input.recal.bam | wc | awk '{ print \$1 }' > $output3.ontarget.txt

          rm -r "$safe_tmp_dir"
        """, "calc_coverage_stats"
     }
    stage_status("calc_coverage_stats", "exit", sample)
}

check_ontarget_perc = {
    stage_status("check_ontarget_perc", "enter", sample)

    var MIN_ONTARGET_PERCENTAGE : 50,
        input_ontarget_file: "qc/${sample}.ontarget.txt" // something about segments messes up $input
    check {
        exec """
            RAW_READ_COUNT=`cat $input_ontarget_file`

            ONTARGET_PERC=`grep -A 1 LIBRARY $input.metrics | tail -1 | awk '{ print int(((\$3 * 2) / "'"$RAW_READ_COUNT"'"))*100 }'`

            [ $RAW_READ_COUNT -eq 0 -o $ONTARGET_PERC -lt $MIN_ONTARGET_PERCENTAGE ]

             """
    } otherwise {
        send text {"On target read percentage for $sample < $MIN_ONTARGET_PERCENTAGE"} to channel: cpipe_operator 
    }
    stage_status("check_ontarget_perc", "exit", sample)
}

calculate_qc_statistics = {
    doc "Calculate additional qc statistics"
    stage_status("calculate_qc_statistics", "enter", sample)
    output.dir="qc"
    // transform("bam") to(sample + ".fragments.tsv") {
    produce("${sample}.fragments.tsv") {
        exec """
            $SAMTOOLS/samtools view $input.recal.bam | python $SCRIPTS/calculate_qc_statistics.py > $output.tsv
        """
    }
    stage_status("calculate_qc_statistics", "exit", sample)
}

gatk_depth_of_coverage = {

    doc "Calculate statistics about depth of coverage for an alignment using GATK"

    stage_status("gatk_depth_of_coverage", "enter", sample)

    output.dir = "qc"
    transform("recal.bam") to(".${target_name}.cov.sample_cumulative_coverage_proportions", 
                         ".${target_name}.cov.sample_interval_statistics") { 
        exec """
            $JAVA -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
               -R $REF
               -T DepthOfCoverage 
               -o $output.sample_cumulative_coverage_proportions.prefix
               --omitDepthOutputAtEachBase
               -I $input.recal.bam
               -ct 1 -ct 10 -ct 20 -ct 50 -ct 100
               -L $target_bed_file.${sample}.bed
        """
    }
    stage_status("gatk_depth_of_coverage", "exit", sample)
}

insert_size_metrics = {

    doc "Generates statistics about distribution of DNA fragment sizes"

    stage_status("insert_size_metrics", "enter", sample)

    var MIN_MEDIAN_INSERT_SIZE : 70,
        MAX_MEDIAN_INSERT_SIZE : 240

    output.dir="qc"
    exec """
        $JAVA -Xmx4g -jar $PICARD_HOME/picard.jar CollectInsertSizeMetrics INPUT=$input.recal.bam O=$output.txt H=$output.pdf
    """

    check {
        exec """
              INSERT_SIZE=`grep -A 1 MEDIAN_INSERT_SIZE $output.txt | cut -f 1 | tail -1 | sed 's/\\.[0-9]*//'`

              echo "Median insert size = $INSERT_SIZE"

              [ $INSERT_SIZE -gt $MIN_MEDIAN_INSERT_SIZE ] && [ $INSERT_SIZE -lt $MAX_MEDIAN_INSERT_SIZE ]
             """, "local"
    } otherwise {
        send text {"""
            WARNING: Insert size distribution for $sample has median out of 
            range $MIN_MEDIAN_INSERT_SIZE - $MAX_MEDIAN_INSERT_SIZE
        """} to channel: cpipe_operator, file: output.pdf
    }
    stage_status("insert_size_metrics", "exit", sample)
}

gap_report = {
    stage_status("gap_report", "enter", sample)

    output.dir="results"
    
    var LOW_COVERAGE_THRESHOLD : 15,
        LOW_COVERAGE_WIDTH : 1,
        input_coverage_file: "qc/${sample}.cov.gz" // something about segments messes up $input

    produce("${run_id}_${sample}.gap.csv") {
        from("$input_coverage_file") {
            exec """
                python $SCRIPTS/gap_annotator.py --min_coverage_ok $LOW_COVERAGE_THRESHOLD --min_gap_width $LOW_COVERAGE_WIDTH --coverage $input_coverage_file --db $BASE/designs/genelists/refgene.txt > $output.csv
            """
        }
    }
    stage_status("gap_report", "exit", sample)
}

summary_report = {
    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    stage_status("summary_report", "enter", sample)

    output.dir="results"

    var input_coverage_file: "qc/${sample}.cov.gz", // something about segments messes up $input
        input_exome_file: "qc/${sample}.exome.gz", 
        input_ontarget_file: "qc/${sample}.ontarget.txt",
        input_fragments_file: "qc/${sample}.fragments.tsv"

    produce("${run_id}_${sample}.summary.htm", "${run_id}_${sample}.summary.md", "${run_id}_${sample}.summary.karyotype.tsv") {
        from("$input_exome_file", "$input_ontarget_file", "$input_fragments_file") {
            exec """
                python $SCRIPTS/qc_report.py --report_cov $input_coverage_file --exome_cov $input_exome_file --ontarget $input_ontarget_file ${inputs.metrics.withFlag("--metrics")} --study $sample --meta $sample_metadata_file --threshold 20 --classes GOOD:95:GREEN,PASS:80:ORANGE,FAIL:0:RED --gc $target_gene_file --gene_cov qc/exon_coverage_stats.txt --write_karyotype $output.tsv --fragments $input_fragments_file --padding $INTERVAL_PADDING_CALL,$INTERVAL_PADDING_INDEL,$INTERVAL_PADDING_SNV > $output.md

                python $SCRIPTS/markdown2.py --extras tables < $output.md | python $SCRIPTS/prettify_markdown.py > $output.htm
            """
        }

        branch.karyotype = output.tsv

        send text {"Sequencing Results for Study $sample"} to channel: cpipe_operator, file: output.htm
    }
    stage_status("summary_report", "exit", sample)
}

summary_report_trio = {
    doc """Generate the summary report and include details of trio members"""
    // this is the quick solution that concatenates the reports together
    stage_status("summary_report_trio", "enter", sample)

    output.dir="results"

    // extract parent sample ids
    def additional_samples = sample_info[sample].pedigree.tokenize(';')[0].tokenize('=')[1].tokenize(',');
    // from_list = add
    def sample_list = additional_samples.collect { "results/${run_id}_${it}.summary.md" }
    sample_list.add("results/${run_id}_${sample}.summary.md")
    stage_status("summary_report_trio", "sample_list", sample_list)
    def sample_string = sample_list.join(' ')
    stage_status("summary_report_trio", "sample_list", sample_string)

    var input_coverage_file: "qc/${sample}.cov.gz", // something about segments messes up $input
        input_exome_file: "qc/${sample}.exome.gz", 
        input_ontarget_file: "qc/${sample}.ontarget.txt",
        input_fragments_file: "qc/${sample}.fragments.tsv"

    produce("${run_id}_${sample}.trio.summary.htm", "${run_id}_${sample}.trio.summary.md") {
        // from("$input_exome_file", "$input_ontarget_file", "$input_fragments_file") {
        from(sample_list) {
            exec """
                cat ${sample_string} > $output.md

                python $SCRIPTS/markdown2.py --extras tables < $output.md | python $SCRIPTS/prettify_markdown.py > $output.htm
            """
        }
    }
 
    stage_status("summary_report_trio", "exit", sample)
}

// this is the final solution that builds a nice report that properly integrates the samples
// currently not complete or used
// summary_report_trio_2 = {
//     doc """Generate the summary report and include details of trio members"""
//     stage_status("summary_report_trio", "enter", sample)
// 
//     output.dir="results"
// 
//     // extract parent sample ids
//     def additional_samples = sample_info[sample].pedigree.tokenize(';')[0].tokenize('=')[1].tokenize(',');
//     // from_list = add
//     def sample_list = additional_samples.collect { "results/${run_id}_${it}.summary.md" }
//     sample_list.add("results/${run_id}_${sample}.summary.md")
//     stage_status("summary_report_trio", "sample_list", sample_list)
//     def sample_string = ' '.join(sample_list)
//     stage_status("summary_report_trio", "sample_list", sample_string)
// 
//     produce("${run_id}_${sample}.trio.summary.htm", "${run_id}_${sample}.trio.summary.md") {
//         from("$input_exome_file", "$input_ontarget_file", "$input_fragments_file") {
//             exec """
//                 python $SCRIPTS/qc_report.py --report_cov $input_coverage_file --exome_cov $input_exome_file --ontarget $input_ontarget_file ${inputs.metrics.withFlag("--metrics")} --study $sample --meta $sample_metadata_file --threshold 20 --classes GOOD:95:GREEN,PASS:80:ORANGE,FAIL:0:RED --gc $target_gene_file --gene_cov qc/exon_coverage_stats.txt --write_karyotype $output.tsv --fragments $input_fragments_file --padding $INTERVAL_PADDING_CALL,$INTERVAL_PADDING_INDEL,$INTERVAL_PADDING_SNV > $output.md
// 
//                 python $SCRIPTS/markdown2.py --extras tables < $output.md | python $SCRIPTS/prettify_markdown.py > $output.htm
//             """
//         }
//     }
//  
//     stage_status("summary_report_trio", "exit", sample)
// }

exon_qc_report = {

    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    stage_status("exon_qc_report", "enter", sample)

    output.dir="results"

    var enable_exon_report : false

    if(!enable_exon_report)  {
        msg "Exon level coverage report not enabled for $target_name"
        return
    }

    produce("${sample}.exon.qc.xlsx", "${sample}.exon.qc.tsv") {
        exec """
             JAVA_OPTS="-Xmx3g" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/exon_qc_report.groovy 
                -cov $input.cov.gz
                -targets $target_bed_file
                -refgene $ANNOVAR_DB/hg19_refGene.txt 
                -x $output.xlsx
                -o $output.tsv
        """
    }
    stage_status("exon_qc_report", "exit", sample)
}

check_coverage = {
    stage_status("check_coverage", "enter", sample)

    output.dir = "qc"

    def medianCov
    transform("cov.gz") to("cov.stats.median", "cov.stats.csv") {

        R {"""
            bam.cov = read.table(pipe("gunzip -c $input.cov.gz"), col.names=c("chr","start","end", "gene", "offset", "cov"))
            meds = aggregate(bam.cov$cov, list(bam.cov$gene), median)
            write.csv(data.frame(Gene=meds[,1],MedianCov=meds$x), "$output.csv", quote=F, row.names=F)
            writeLines(as.character(median(bam.cov$cov)), "$output.median")
        """}

        // HACK to ensure file sync on distributed file system
        file(output.dir).listFiles()
        medianCov = Math.round(file(output.median).text.toFloat())
    }

    check {
        exec "[ $medianCov -ge $MEDIAN_COVERAGE_THRESHOLD ]"
    } otherwise {
        // It may seem odd to call this a success, but what we mean by it is that
        // Bpipe should not fail the whole pipeline, merely this branch of it
        succeed report('templates/sample_failure.html') to channel: cpipe_operator, 
                                                           median: medianCov, 
                                                           file:output.csv, 
                                                           subject:"Sample $sample has failed with insufficient median coverage ($medianCov)"
    }
    stage_status("check_coverage", "exit", sample)
}

check_karyotype = {

    doc "Compare the inferred sex of the sample to the inferred karyotype from the sequencing data"
    stage_status("check_karyotype", "enter", sample)

    def karyotype_file = "results/" + run_id + '_' + sample + '.summary.karyotype.tsv'
    check {
        exec """
            [ `grep '^Sex' $karyotype_file | cut -f 2` == "UNKNOWN" ] || [ `grep '^Sex' $karyotype_file | cut -f 2` == `grep 'Inferred Sex' $karyotype_file | cut -f 2` ]
        """
    } otherwise {
        // It may seem odd to call this a success, but what we mean by it is that
        // Bpipe should not fail the whole pipeline, merely this branch of it
        //succeed report('templates/sample_failure.html') to channel: cpipe_operator, 
        //                                                   median: medianCov, 
        //                                                   file: karyotype_file,
        //                                                   subject:"Sample $sample has a different sex than inferred from sequencing data"
        stage_status("check_karyotype", "gender comparison failed", sample)
    }
    stage_status("check_karyotype", "exit", sample)
}

qc_excel_report = {

    doc "Create an excel file containing a summary of QC data for all the samples for a given target region"
    stage_status("qc_excel_report", "enter", sample)

    var LOW_COVERAGE_THRESHOLD : 15,
        LOW_COVERAGE_WIDTH : 1

    output.dir="results"

    def samples = sample_info.grep { it.value.target == target_name }.collect { it.value.sample }
    produce(target_name + ".qc.xlsx") {
            exec """
                JAVA_OPTS="-Xmx16g -Djava.awt.headless=true" $GROOVY -cp $GROOVY_NGS/groovy-ngs-utils.jar:$EXCEL/excel.jar $SCRIPTS/qc_excel_report.groovy 
                    -s ${target_samples.join(",")} 
                    -t $LOW_COVERAGE_THRESHOLD
                    -w $LOW_COVERAGE_WIDTH
                    -low qc ${inputs.dedup.metrics.withFlag('-metrics')}
                    -o $output.xlsx
                    -p $run_id
                    $inputs.sample_cumulative_coverage_proportions  
                    $inputs.sample_interval_statistics 
                    $inputs.gz
            ""","qc_excel_report"
    }
    stage_status("qc_excel_report", "exit", sample)
}

provenance_report = {
    stage_status("provenance_report", "enter", sample)
    branch.sample = branch.name
    output.dir = "results"
    produce(run_id + '_' + sample + ".provenance.pdf") {
       send report("scripts/provenance_report.groovy") to file: output.pdf
    }
    stage_status("provenance_report", "exit", sample)
}

filtered_on_exons = {
    doc "Create a bam filtered on exons with 100bp padding and excluding the incidentalome"
    stage_status("filtered_on_exons", "enter", sample)
    // bedtools exons.bed + padding100bp - incidentalome
    // TODO this might be faster if we sorted the bam and used -sorted
    var GENE_BAM_PADDING: 100

    def safe_tmp = ['tmp', UUID.randomUUID().toString()].join( '' )

    output.dir = "results"

    produce("${run_id}_${branch.name}.filtered_on_exons.bam") {
        exec """
            python $SCRIPTS/filter_bed.py --include $BASE/designs/genelists/incidentalome.genes.txt < $BASE/designs/genelists/exons.bed |
            $BEDTOOLS/bin/bedtools slop -g $HG19_CHROM_INFO -b $GENE_BAM_PADDING -i - > $safe_tmp 
            
            python $SCRIPTS/filter_bed.py --exclude $BASE/designs/genelists/incidentalome.genes.txt < $BASE/designs/genelists/exons.bed |
            $BEDTOOLS/bin/bedtools slop -g $HG19_CHROM_INFO -b $GENE_BAM_PADDING -i - | 
            $BEDTOOLS/bin/bedtools subtract -a - -b $safe_tmp | 
            sort -k1,1 -k2,2n |
            $BEDTOOLS/bin/bedtools intersect -a $input.recal.bam -b stdin > $output.bam
    
            rm "$safe_tmp"
        """
    }
    stage_status("filtered_on_exons", "exit", sample)
}

variant_filtering_report = {
    doc "generate a report of all variants and where they were filtered"
    output.dir = "variants"
    produce("variant_filtering_report.tsv") {
        exec """
            python $SCRIPTS/variant_filtering.py  --source_dir $output.dir > $output
        """
    }
}

// no longer part of the default pipeline
variant_bams = {
    doc "Create a bam file for each variant containing only reads overlapping 100bp either side of that variant"
    stage_status("variant_bams", "enter", sample)

    output.dir = "results/variant_bams"

    from("${run_id}_${branch.name}*.lovd.tsv", "${branch.name}.*.recal.bam") {   
        // Slight hack here. Produce a log file that bpipe can track to confirm that the bams were produced.
        // Bpipe is not actually tracking the variant bams themselves. 
        // TODO should these be generated from the HC output?
        produce(branch.name + ".variant_bams_log.txt") {
            exec """
                python $SCRIPTS/variant_bams.py --bam $input.recal.bam --tsv $input.tsv --outdir $output.dir --log $output.txt --samtoolsdir $SAMTOOLS
            """
        }
    }
    stage_status("variant_bams", "exit", sample)
}

///////////////////////////////////////////////////////////////////
// segments
///////////////////////////////////////////////////////////////////
analysis_ready_reports = segment {
    // parallel doesn't work properly here
    calc_coverage_stats + 
    check_ontarget_perc + 
    calculate_qc_statistics + 
    summary_report + 
    exon_qc_report + 
    gap_report +
    gatk_depth_of_coverage +
    insert_size_metrics +
    filtered_on_exons + index_bam 
}

analysis_ready_checks = segment {
    check_coverage +
    check_karyotype
}

