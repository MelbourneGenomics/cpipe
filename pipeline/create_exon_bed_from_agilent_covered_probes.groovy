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
load '../config.groovy'

@produce("CARDIAC101.genes.bed")
covered_probes_to_genes = {
    output.dir = "work"
    exec """
	agilent_covered_probes_to_genes  $input.bed > $output.bed
    """
}

@produce("CARDIAC101.exons.bed")
create_exon_bed = {
    output.dir = "work"
    exec """
	create_exon_bed -c $input.bed $BASE/tools/annovar/humandb/hg19_refGene.txt $input.txt $output.bed
    """
}

@filter("sort")
sort_bed = {
    output.dir = "work"
    exec """
	bedtools sort -i $input.bed > $output.bed
    """
}


inputs "bed" : "Agilent design file of covered regions",
       "txt" : "VCGS transcripts file containing list of prioritized transcripts"

run {
    covered_probes_to_genes + create_exon_bed + sort_bed
}
