#!/usr/bin/env python3
'''
###########################################################################
#
# This file is part of Cpipe.
#
# Cpipe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, under version 3 of the License, subject
# to additional terms compatible with the GNU General Public License version 3,
# specified in the LICENSE file that is part of the Cpipe distribution.
#
# Cpipe is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.  If not, see <http:#www.gnu.org/licenses/>.
#
###########################################################################

#
# Reads a file of exons (arg1) and the genome annotation file from annovar
# (arg2) and outputs all the variants that are within <n> base pairs of an exon
# (arg3)
'''

import csv
import sys

def process(exons, genome, width):
    '''
      write variants to stdout
    '''
    search_width = int(width)
    exon_file = csv.reader(open(exons), delimiter='\t')

    exons = {}
    count = 0
    for ex in exon_file:
        chrom, start, end, desc = (ex[0], ex[1], ex[2], ex[3])
        exons.setdefault(chrom, []).append([int(start), int(end), desc])
        count += 1

    #print "Parsed %d exons (%d chromosomes)" % (count, len(exons))

    wout = csv.writer(sys.stdout)

    # Now read the annovar file
    annovar_file = csv.reader(open(genome))
    #header = annovar_file.next()

    #w.writerow(header)

    for line in annovar_file:
        #gene = l[1]
        #if gene != 'CACNB2':
        #     continue

        chrom = line[21]
        start = int(line[22])
        end = int(line[23])
        #aachange = line[3]

        #print  ",".join([chrom, str(start), str(end), aachange])

        if chrom not in exons:
            print("WARNING: Chromosome not in capture in output variants: " + chrom, file=sys.stderr)
            continue

        for exon_start, exon_end, desc in exons[chrom]:
            #if exon_start == 18439811 and start == 18439809:
            #print "Exon [%d,%d] " % (exon_start, exon_end)
            if end <= exon_start and end > exon_start - search_width:
                line[0] = "extra_splicing;"+line[0]
                wout.writerow(line)
            elif start >= exon_end and start < exon_end + search_width:
                line[0] = "extra_splicing;"+line[0]
                wout.writerow(line)
            else:
                # Not interesting
                pass

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python add_splice_variants.py <exons.bed> <genome_summary.csv> <width>")
        sys.exit(1)

    process(sys.argv[1], sys.argv[2], sys.argv[3])
