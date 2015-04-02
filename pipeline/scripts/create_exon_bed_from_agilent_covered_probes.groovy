load '../config.groovy'

@produce("CARDIAC101.genes.bed")
covered_probes_to_genes = {
    output.dir = "work"
    exec """
	groovy $SCRIPTS/agilent_covered_probes_to_genes  $input.bed > $output.bed
    """
}

@produce("CARDIAC101.exons.bed")
create_exon_bed = {
    output.dir = "work"
    exec """
	python $SCRIPTS/create_exon_bed.py -c $input.bed $BASE/tools/annovar/humandb/hg19_refGene.txt $input.txt $output.bed
    """
}

@filter("sort")
sort_bed = {
    output.dir = "work"
    exec """
	$BEDTOOLS/bin/bedtools sort -i $input.bed > $output.bed
    """
}


inputs "bed" : "Agilent design file of covered regions",
       "txt" : "VCGS transcripts file containing list of prioritized transcripts"

run {
    covered_probes_to_genes + create_exon_bed + sort_bed
}
