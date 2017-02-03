#!/usr/bin/env python3
# vim: sw=4:expandtab:cindent:ts=4
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
# Reads the output summary from Annovar and adds a column to indicate 
# whether the mutation in the row affects one of the transcripts that
# are identified in the prioritised transcripts file.
#
###########################################################################

import csv
import sys
import re

def main():
    # Transcripts of interest
    transcripts = sys.argv[1]

    # Summary exome_summary.csv file from Annovar
    summary = sys.argv[2]

    # Full exonic_variant_function file from Annovar
    full = sys.argv[3]

    # Read all the transcripts of interest
    txs = {}
    for i in open(sys.argv[1]):
        txs[re.sub('\\.[0-9].*$', '', i.strip())] = True

    w = csv.writer(sys.stdout)

    # Read the CSV summary file and examine each variant
    reader = csv.reader(open(summary))

    # First read the header and since by default some columns don't have headers,
    # as a side benefit we fix those
    header = next(reader)

    # Fix missing column headings
    header.append('PRIORITY_TX')

    chr_index = header.index('Chr')
    pos_index = header.index('Start')
    gene_index = header.index('Gene')
    aachange_index = header.index('AAChange')

    w.writerow(header)

    debug = False

    for l in reader:
        gene = l[gene_index]
        chr = l[chr_index]
        start = l[pos_index]
        aachange = l[aachange_index]

        if gene == 'unknown':
            continue

        # Search for the aachange in the full file
        # to get the full list of isoforms / transcripts and see if any are
        # flagged as of interest
        f = open(full)
        row = l
        found_vtx = ''
        for v in csv.reader(f, delimiter='\t'):
            # print >>sys.stderr, str(v)

            if v[1] == 'unknown':
                continue

            # print >>sys.stderr, "Check: %s:%s vs %s:%s" % (chr, start, v[3],v[4])

            # Same location
            if v[3] != chr:
                continue

            if v[4] != start:
                continue

            # Has to be the same gene
            if gene not in [x.split(':')[0] for x in v[2].split(',')]:
                continue

            # Column 2 is in the following format:
            # TTN:NM_003319:exon73:c.G18133A:p.D6045N,TTN:NM_133432:exon74:c.G18508A:p.D6170N
            # We want to report only the transcript and the AA change

            try:
                vtxs = list(zip([x.split(':')[1] for x in v[2].split(',')],  # the NM_.. transcript id
                                [':'.join(x.split(':')[3:5]) for x in v[2].split(',')]))  # the AA change

                if debug:
                    print("Full transcripts for %s are %s" % (aachange, vtxs), file=sys.stderr)

                # vtxs is a list of transcripts, each element is another list
                # of 2 elements, (tx name, aa change)
                vtxs_flag = [x for x in vtxs if x[0] in txs]
                if vtxs_flag:
                    if debug:
                        print("Full transcripts for %s are %s" % (aachange, vtxs))
                        print("FLAGGED: %s" % vtxs_flag)
                    found_vtx = ";".join([":".join(f) for f in vtxs_flag])
            except Exception as e:
                found_vtx = 'Error: Please check manually'

        row.append(found_vtx)
        w.writerow(row)
        f.close()


if __name__ == '__main__':
    main()
