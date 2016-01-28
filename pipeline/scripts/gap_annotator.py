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
# Purpose:
#   Generate annotated gap CSV
# Usage:
#   See main() for arguments
#
###########################################################################
'''

import datetime
import sys

def write_log(log, msg):
    '''
        write a date stamped message to log
    '''
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    if log is not None:
        log.write('%s: %s\n' % (now, msg))

def median(items):
    '''
      return the median of a list of items
    '''
    sorted_list = sorted(items)
    if len(items) % 2 == 0:
        high = len(items) / 2
        return (sorted_list[high] + sorted_list[high-1]) / 2.
    else:
        mid = (len(items) - 1) / 2
        return sorted_list[mid]

def annotate_gap(gap, target):
    '''
      write out the found gap
    '''
    # note that start and end are inclusive
    target.write('{0},{1},{2},{3},{4},{5},{6},{7}\n'.format(gap['chr'], gap['gene'], gap['start'] + gap['start_offset'] - 1, gap['start'] + gap['start_offset'] + gap['length'] - 1 - 1, min(gap['coverage']), max(gap['coverage']), median(gap['coverage']), gap['length']))

def find_gaps(coverage, min_width, max_coverage, target, log):
    '''
        find gaps and annotate
    '''
    target.write('{0},{1},{2},{3},{4},{5},{6},{7}\n'.format('Chr', 'Gene', 'Start', 'End', 'Min Cov', 'Max Cov', 'Median Cov', 'Width'))
    current = None
    gaps = 0
    idx = 0
    write_log(log, 'finding gaps...')
    for idx, line in enumerate(coverage):
        fields = line.strip('\n').split('\t') # chr, start, end, gene, offset, coverage
        if len(fields) > 5:
            coverage = int(fields[5])
            if current is None: # not in gap
                if coverage <= max_coverage: # start a new gap
                    current = {'start': int(fields[1]), 'start_offset': int(fields[4]), 'length': 1, 'chr': fields[0], 'gene': fields[3], 'coverage': [coverage]}
            else: # in gap
                if fields[0] == current['chr'] and int(fields[1]) == current['start'] and fields[3] == current['gene'] and coverage <= max_coverage: # continue gap
                    current['length'] += 1
                    current['coverage'].append(coverage)
                else: # end of gap
                    if current['length'] >= min_width: # write it out
                        annotate_gap(current, target)
                        gaps += 1
                    if coverage <= max_coverage: # start a new gap?
                        current = {'start': int(fields[1]), 'start_offset': int(fields[4]), 'length': 1, 'chr': fields[0], 'gene': fields[3], 'coverage': [coverage]}
                    else:
                        current = None
        else:
            write_log(log, 'skipped line {0}'.format(idx))

        if idx % 100000 == 0:
            write_log(log, 'finding gaps: {0} lines processed; {1} gaps found...'.format(idx, gaps))

    # still a gap in progress?
    if current is not None and current['length'] >= min_width:
        annotate_gap(current, target)
    write_log(log, 'finding gaps: {0} lines; {1} gaps: done'.format(idx, gaps))

def main():
    '''
        parse command line and execute
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Generate gap report')
    parser.add_argument('--min_coverage_ok', required=False, type=int, default=0, help='maximum value to consider to be low coverage')
    parser.add_argument('--min_gap_width', required=False, type=int, default=1, help='minimum width of a gap to report')
    parser.add_argument('--coverage', required=True, help='coverage file to examine for gaps')
    args = parser.parse_args()
    find_gaps(open(args.coverage, 'r'), args.min_gap_width, args.min_coverage_ok, sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()
