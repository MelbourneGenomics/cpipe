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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.    If not, see <http:#www.gnu.org/licenses/>.
#
###########################################################################
# Purpose:
#   filter tab separated output based on the transcript:
#   remove XMs that have an NM at the same location
# Usage:
#   filter_transcripts.py < variants.tsv 1>filtered.tsv
####################################################################################
'''

import argparse
import collections
import datetime
import sys

# which column name contains the transcript
TRANSCRIPT_COLUMN = 'Feature'

# how we identify a variant
IDENTIFYING_COLUMNS = ['CHROM', 'POS', 'REF', 'ALT']

def filter_tsv(instream, outstream, log):
    '''
        read instream and filter out XM transcripts that are in a location that 
        already has an NM transcript.
        write filtered transcripts to outstream.
    '''
    # get header
    first = True
    include = collections.defaultdict(list)
    has_nm = set()
    log.write('reading...\n')
    count = 0
    counts = {'nm_orig': 0, 'xm_orig': 0, 'other_orig': 0, 'nm_new': 0, 'xm_new': 0}
    for count, line in enumerate(instream):
        if first: 
            # get header info
            first = False
            header_line = line
            header = line.strip('\n').split('\t')
            if TRANSCRIPT_COLUMN not in header:
                log.write('ERROR: {0} not found in header\n'.format(TRANSCRIPT_COLUMN))
                return
            transcript_idx = header.index(TRANSCRIPT_COLUMN)
            identifiers = []
            for identifier in IDENTIFYING_COLUMNS:
                if identifier not in header:
                    log.write('ERROR: {0} not found in header\n'.format(identifier))
                    return
                identifiers.append(header.index(identifier))
        else:
            fields = line.strip('\n').split('\t')
            # build identifier
            identifier = '\t'.join([fields[x] for x in identifiers])
            transcript = fields[transcript_idx]
            if transcript.startswith('NM_'):
                # if the transcript starts with NM and NM hasn't been seen before,
                # all previously seen XMs are removed from the list
                if identifier not in has_nm:
                    new_list = [ line ]
                    for old_item in include[identifier]:
                        old_item_transcript = old_item.strip('\n').split('\t')[transcript_idx]
                        if not old_item_transcript.startswith('XM_'):
                            new_list.append(old_item)
                    include[identifier] = new_list
                # if the transcript starts with NM and we've already seen an NM, 
                # add the NM
                else: 
                    include[identifier].append(line)
                has_nm.add(identifier) 
                counts['nm_orig'] += 1
            elif transcript.startswith('XM_'):
                # if the transcript starts with XM only add it if we haven't
                # seen an nm
                if not identifier in has_nm:
                    include[identifier].append(line)
                counts['xm_orig'] += 1
            # anything that's not an XM or NM is added
            else: 
                include[identifier].append(line)
                counts['other_orig'] += 1
    log.write('read {0} lines\n'.format(count))

    # write out filtered list
    count = 0
    outstream.write(header_line)
    for identifier in include:
        for line in include[identifier]:
            transcript = line.strip('\n').split('\t')[transcript_idx]
            if transcript.startswith('NM_'):
                counts['nm_new'] += 1
            elif transcript.startswith('XM_'):
                counts['xm_new'] += 1
            outstream.write(line)
            count += 1
    log.write('done writing {0} lines. stats: {1}\n'.format(count, counts))

def main():
    '''
        run from command line
    '''
    parser = argparse.ArgumentParser(description='Filter TSV')
    args = parser.parse_args()
    filter_tsv(sys.stdin, sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()
