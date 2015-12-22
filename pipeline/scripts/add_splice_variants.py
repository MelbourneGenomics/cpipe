#!/usr/bin/env python
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

import csv, sys

if len(sys.argv) < 4:
    print "Usage: python add_splice_variants.py <exons.bed> <genome_summary.csv> <width>"
    sys.exit(1)

search_width = int(sys.argv[3])

exon_file = csv.reader(open(sys.argv[1]), delimiter='\t')

exons = {}
count = 0
for e in exon_file:
    chr, start, end, desc = (e[0], e[1], e[2], e[3])
    exons.setdefault(chr, []).append([int(start), int(end), desc])
    count += 1

#print "Parsed %d exons (%d chromosomes)" % (count, len(exons))

w = csv.writer(sys.stdout)

# Now read the annovar file
annovar_file = csv.reader(open(sys.argv[2]))
header = annovar_file.next()

#w.writerow(header)

for l in annovar_file:
    gene = l[1]
    #if gene != 'CACNB2':
    #     continue

    chr = l[21]
    start = int(l[22])
    end = int(l[23])
    aachange = l[3]

    #print  ",".join([chr, str(start), str(end), aachange])

    if not chr in exons:
        print >>sys.stderr, "WARNING: Chromosome not in capture in output variants: " + chr
        continue

    for exon_start, exon_end, desc in exons[chr]:
        #if exon_start == 18439811 and start == 18439809:
        #print "Exon [%d,%d] " % (exon_start, exon_end)
        if end <= exon_start and end > exon_start - search_width:
            l[0] = "extra_splicing;"+l[0]
            w.writerow(l)
        elif start >= exon_end and start < exon_end + search_width:
            l[0] = "extra_splicing;"+l[0]
            w.writerow(l)
        else:
            # Not interesting
            pass

