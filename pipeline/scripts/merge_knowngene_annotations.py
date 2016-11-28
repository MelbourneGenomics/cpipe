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
# Reads output summaries from Annovar based on RefSeq and UCSC respectively and for 
# each row in one, tries to find the same variant in the other.  A new column is
# added showing the alternative annotation next to the original
# where such an annotation can be located.
#
# This was added to mitigate against problems of incorrect annotations in RefSeq.
# In general it was found that in the rare cases when RefSeq was incorrect, UCSC
# was correct (or at least, inconsistent).  Therefore a way to at least flag 
# potential for problems was to annotate with both, and flag when there is a 
# difference between the two.
#
import csv,sys,re


# Summary exome_summary.csv file from Annovar
summary = sys.argv[1]

# The alternative annotations
alt = sys.argv[2]

w = csv.writer(sys.stdout)

# Read the CSV summary file and examine each variant
reader = csv.reader(open(summary))

# First read the header and since by default some columns don't have headers,
# as a side benefit we fix those
header = next(reader)

VCGS_TX_COL = header.index("VCGS_TX")

header = header[0:3] + ["AAChange_RefSeq", "AAChange_UCSC"] + header[4:] + ["AA_Match"]

w.writerow(header)

debug = False

def find_alt_annotation(gene,chr,start,obs, aa):
        altf = open(alt)
        altr = csv.reader(altf)
        try:
            for l in altr:
                    if gene == l[1] and chr == l[21] and start == l[22] and obs == l[25]:
                            # Found the variant (we hope)
                            if debug:
                                    print("Matched %s:%s:%s:%s with aa change %s vs %s" % (gene,chr,start,obs, aa, l[VCGS_TX_COL]), file=sys.stderr)

                            return l[3]
            return ""
        finally:
                altf.close()
                            


for l in reader:
        gene = l[1]
        chr = l[21]
        start = l[22]
        aachange = l[3]

        if gene == 'unknown':
                continue

        # Search for the location in the first file
        alt_annotation = find_alt_annotation(gene,chr,start,l[25], l[VCGS_TX_COL])

        change1 = re.sub('^.*?:', '', l[3]) 
        change2 = re.sub('^.*?:', '', alt_annotation) 

        # The above are in the form c.G14744A:p.R4915H
        # However we want to ignore positional changes and only focus
        # on amino acid inconsistencies.  This is because we don't have a good 
        # way to translate a RefSeq protein identifier (NM_...) to a UCSC one 
        # (uc.xxxxx).  Hence we may be comparing different isoforms.
        #
        # Hacky way to do that is to remove all the numbers and then
        # compare strings

        aa_match = re.sub('[0-9]*', '', change1) == re.sub('[0-9]*', '', change2)

        l = l[0:4] + [alt_annotation] + l[4:] + [aa_match]

        w.writerow(l)

