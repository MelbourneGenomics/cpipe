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
# Filter a coverage file for regions < a threshold in coverage and
# output a unique identifier as an annotation for those regions
# This eases further processing in R
###########################################################################

import sys, csv

if len(sys.argv) < 4:
    print "\nUsage: %s <coverage file> <output file> <coverage threshold>\n" % sys.argv[0]
    exit(1)

cov = csv.reader(open(sys.argv[1]), delimiter='\t')
output = csv.writer(open(sys.argv[2],'wb'), delimiter='\t')

coverage_threshold = int(sys.argv[3])

in_block = False
block_index = 0
block_length = -1
low_count = 0
high_count = 0
total_count = 0
prev_info = None

for line in cov:
    chr,exon_start,exon_end,gene = line[0:4]
    info = gene + exon_start
    offset,cov = line[-2::]
    exon_start,exon_end,offset,cov = int(exon_start),int(exon_end),int(offset),int(cov)

    if prev_info != info:
        in_block = False

    prev_info = info

    total_count += 1

    if cov < coverage_threshold:
        #print "offset = %d cov = %d" % (offset,cov)
        if not in_block:
            in_block = True
            block_index += 1
            block_length = 1

        line.append(block_index)
        line.append(block_length)
        output.writerow(line)
        low_count += 1
        block_length += 1
    else:
        in_block = False
        high_count += 1

print "Found %d areas with coverage < %d (total %d bases above threshold + %d below = %d total)" % (block_index, coverage_threshold, high_count, low_count, total_count)
