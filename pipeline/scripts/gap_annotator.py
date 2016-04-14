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

import collections
import datetime
import math
import os
import sys

# interval matching from https://bitbucket.org/james_taylor/bx-python/raw/ebf9a4b352d3/lib/bx/intervals/operations/quicksect.py
import random

class IntervalTree(object):
    '''
         fast interval finder
    '''
    def __init__(self):
        '''
            dictionary of intervals
        '''
        self.chroms = {}

    def insert(self, interval, linenum=0, other=None):
        '''
            add a new interval
        '''
        chrom = interval.chrom
        start = interval.start
        end = interval.end
        if interval.chrom in self.chroms:
            self.chroms[chrom] = self.chroms[chrom].insert(start, end, linenum, other)
        else:
            self.chroms[chrom] = IntervalNode(start, end, linenum, other)

    def intersect(self, interval, report_func):
        '''
            add a new interval
        '''
        chrom = interval.chrom
        start = interval.start
        end = interval.end
        if chrom in self.chroms:
            self.chroms[chrom].intersect(start, end, report_func)

    def traverse(self, func):
        '''
            iterate over all intervals
        '''
        for item in self.chroms.itervalues():
            item.traverse(func)

class IntervalNode(object):
    '''
        an interval in a tree
    '''

    def __init__(self, start, end, linenum=0, other=None):
        '''
        Python lacks the binomial distribution, so we convert a
        uniform into a binomial because it naturally scales with
        tree size.  Also, python's uniform is perfect since the
        upper limit is not inclusive, which gives us undefined here.
        '''
        self.priority = math.ceil((-1.0 / math.log(.5))* math.log(-1.0 / (random.uniform(0, 1)- 1)))
        self.start = start
        self.end = end
        self.maxend = self.end
        self.minend = self.end
        self.left = None
        self.right = None
        self.linenum = linenum
        self.other = other

    def insert(self, start, end, linenum=0, other=None):
        '''
            add an interval to the tree
        '''
        root = self
        if start > self.start:
            # insert to right tree
            if self.right:
                self.right = self.right.insert(start, end, linenum, other)
            else:
                self.right = IntervalNode(start, end, linenum, other)
            # rebalance tree
            if self.priority < self.right.priority:
                root = self.rotateleft()
        else:
            # insert to left tree
            if self.left:
                self.left = self.left.insert(start, end, linenum, other)
            else:
                self.left = IntervalNode(start, end, linenum, other)
            # rebalance tree
            if self.priority < self.left.priority:
                root = self.rotateright()
        if root.right and root.left:
            root.maxend = max(root.end, root.right.maxend, root.left.maxend)
            root.minend = min(root.end, root.right.minend, root.left.minend)
        elif root.right:
            root.maxend = max(root.end, root.right.maxend)
            root.minend = min(root.end, root.right.minend)
        elif root.left:
            root.maxend = max(root.end, root.left.maxend)
            root.minend = min(root.end, root.left.minend)
        return root

    def rotateright(self):
        '''
            rearrange tree
        '''
        root = self.left
        self.left = self.left.right
        root.right = self
        if self.right and self.left:
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend)
        return root

    def rotateleft(self):
        '''
            rearrange tree
        '''
        root = self.right
        self.right = self.right.left
        root.left = self
        if self.right and self.left:
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend)
        return root

    def intersect(self, start, end, report_func):
        '''
            find intersection
        '''
        if start < self.end and end > self.start:
            report_func(self)
        if self.left and start < self.left.maxend:
            self.left.intersect(start, end, report_func)
        if self.right and end > self.start:
            self.right.intersect(start, end, report_func)

    def traverse(self, func):
        '''
            find all
        '''
        if self.left:
            self.left.traverse(func)
        func(self)
        if self.right:
            self.right.traverse(func)

def traversal_handler_builder(start, end):
    '''
        helper to traverse all
    '''
    data = [None, 3e9]
    def traversal_handler(interval):
        '''
            closure for traversal
        '''
        if data[0] is None:
            data[0] = interval
            if interval.start > end:
                data[1] = interval.start - end
            else:
                data[1] = start - interval.end
        elif interval.start > end and interval.start - end < data[1]:
            data[0] = interval
            data[1] = interval.start - end
        elif start > interval.end and start - interval.end < data[1]:
            data[0] = interval
            data[1] = start - interval.end
    return data, traversal_handler

def annotate_gap(gap, data_source, log):
    '''
        given a gap, annotate it
    '''
    # find overlapping cds interval, or nearest
    target = []
    if gap['chr'] in data_source['cds'].chroms:
        start_pos = gap['start'] + gap['start_offset'] - 1
        end = start_pos + gap['length']
        data_source['cds'].chroms[gap['chr']].intersect(start_pos, end, lambda x: target.append(x))# find overlap
        if len(target) == 0: # no overlap
            # find closest by traversing the entire tree
            result, handler = traversal_handler_builder(start_pos, end)
            data_source['cds'].chroms[gap['chr']].traverse(handler)
            return {'interval': result[0], 'distance': result[1]}
        else: # overlap found
            if len(target) > 1:
                #write_log(log, "annotate_gap: {0} overlaps found for {1}".format(len(target), gap['start']))
                pass
            #write_log(log, "annotate_gap: {0} {1} {2}".format(target[0].start, target[0].end, target[0].end - target[0].start))# debug

            # intersect the gap with the coding region
            #for i in xrange(0, len(target)):
            #    intersect = find_intersect(target[i].start, target[i].end, start_pos, end)
            #    write_log(log, '{0}. intersect {1} gap {2}->{3} cds {4}->{5}'.format(i, intersect, start_pos, end, target[i].start, target[i].end))

            # what is the intersection relative to the coding region?
            intersect = find_intersect(target[0].start, target[0].end, start_pos, end)
            coding_intersect = [intersect[0] - target[0].start + 1, intersect[1] - target[0].start] # start=1 based, end=inclusive
            #if target[0].other['strand'] == '+': # count from start
            #    coding_intersect = [x - target[0].start for x in intersect]
            #else: # count from end
            #    coding_intersect = [target[0].end - x for x in intersect]
            codon_positions = [int(math.ceil(x / 3.0)) for x in coding_intersect]

            # determine exon rank
            if target[0].other['strand'] == '+':
                exon_rank = target[0].other['number']
            else:
                exon_rank = target[0].other['count'] - target[0].other['number'] + 1

            #write_log(log, "annotate_gap: coding region intersect: {0} codons: {1}".format(coding_intersect, codon_positions))
            return {'interval': target[0], 'distance': 0, 'coding_intersect': coding_intersect, 'codon_positions': codon_positions, 'rank': exon_rank}
    else:
        return None # unexpected chromosome

def write_log(log, msg):
    '''
        write a date stamped message to log
    '''
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    if log is not None:
        log.write('%s: %s\n' % (now, msg))

def run(cmd, log):
    '''
        executes a system command
    '''
    write_log(log, 'executing {0}\n'.format(cmd))
    os.system(cmd)

def median(items):
    '''
      return the median of a list of items
    '''
    sorted_list = sorted(items)
    if len(items)% 2 == 0:
        high = len(items)/ 2
        return (sorted_list[high] + sorted_list[high-1])/ 2.
    else:
        mid = (len(items)- 1)/ 2
        return sorted_list[mid]

def write_gap(gap, target, data_source, log):
    '''
      write out the found gap
    '''
    # note that start and end are inclusive
    annotation = annotate_gap(gap, data_source, log)
    if annotation is None: # shouldn't happen unless things are really wrong
        target.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16}\n'.format(gap['chr'], gap['gene'], gap['start'] + gap['start_offset'] - 1, gap['start'] + gap['start_offset'] + gap['length'] - 1 - 1, min(gap['coverage']), max(gap['coverage']), median(gap['coverage']), gap['length'], 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A'))
    elif annotation['distance'] == 0: # gap overlapping coding sequence
        target.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16}\n'.format(gap['chr'], gap['gene'], gap['start'] + gap['start_offset'] - 1, gap['start'] + gap['start_offset'] + gap['length'] - 1 - 1, min(gap['coverage']), max(gap['coverage']), median(gap['coverage']), gap['length'], annotation['interval'].other['name'], annotation['interval'].other['strand'], annotation['distance'], annotation['coding_intersect'][0], annotation['coding_intersect'][1], annotation['codon_positions'][0], annotation['codon_positions'][1], annotation['interval'].other['number'], annotation['rank']))
    else: # nearest distance
        target.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16}\n'.format(gap['chr'], gap['gene'], gap['start'] + gap['start_offset'] - 1, gap['start'] + gap['start_offset'] + gap['length'] - 1 - 1, min(gap['coverage']), max(gap['coverage']), median(gap['coverage']), gap['length'], annotation['interval'].other['name'], annotation['interval'].other['strand'], annotation['distance'], 'N/A', 'N/A', 'N/A', 'N/A', annotation['interval'].other['number'], 'N/A'))

def find_gaps(coverage, min_width, max_coverage, target, data_source, log):
    '''
        find gaps and annotate
    '''
    target.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16}\n'.format('Chr', 'Gene', 'Start', 'End', 'Min Cov', 'Max Cov', 'Median Cov', 'Width', 'Tx Name', 'Strand', 'CDS Distance', 'CDS Overlap Start', 'CDS Overlap End', 'AA Overlap Start', 'AA Overlap End', 'Exon Number', 'Exon Rank'))
    current = None
    gaps = 0
    idx = 0
    write_log(log, 'finding gaps...')
    for idx, line in enumerate(coverage):
        fields = line.strip('\n').split('\t')# tab separated: chr, start, end, gene, offset, coverage
        if len(fields) > 5 and fields[0].lower() != 'chr':
            coverage = int(fields[5])
            if current is None: # not in gap
                if max_coverage == -1 or coverage <= max_coverage: # start a new gap
                    current = {'start': int(fields[1]), 'start_offset': int(fields[4]), 'length': 1, 'chr': fields[0], 'gene': fields[3], 'coverage': [coverage]}
            else: # in gap
                if fields[0] == current['chr'] and int(fields[1]) == current['start'] and fields[3] == current['gene'] and (max_coverage == -1 or coverage <= max_coverage): # continue gap
                    current['length'] += 1
                    current['coverage'].append(coverage)
                else: # end of gap
                    if current['length'] >= min_width: # write it out
                        write_gap(current, target, data_source, log)
                        gaps += 1
                    if max_coverage == -1 or coverage <= max_coverage: # start a new gap?
                        current = {'start': int(fields[1]), 'start_offset': int(fields[4]), 'length': 1, 'chr': fields[0], 'gene': fields[3], 'coverage': [coverage]}
                    else:
                        current = None
        else:
            write_log(log, 'skipped line {0}'.format(idx))

        if idx % 100000 == 0:
            write_log(log, 'finding gaps: {0} lines processed; {1} gaps found...'.format(idx, gaps))

    # still a gap in progress?
    if current is not None and current['length'] >= min_width:
        write_gap(current, target, data_source, log)
        gaps += 1
    write_log(log, 'finding gaps: {0} lines; {1} gaps: done'.format(idx, gaps))

Interval = collections.namedtuple('Interval', ['start', 'end', 'chrom'])

def find_intersect(candidate_start, candidate_end, range_start, range_end):
    '''
        find the overlapping region
    '''
    if candidate_end <= range_start or candidate_start >= range_end:
        return None
    else:
        return (max(candidate_start, range_start), min(candidate_end, range_end))

def init_db(target, log):
    '''
        prepare annotation db
    '''
    write_log(log, 'starting init_db...')
    result = {'cds': IntervalTree()}
    added = 0
    first = True
    for i, line in enumerate(target):
        if first:
            first = False
            continue
        fields = line.strip('\n').split('\t')
        chrom = fields[2]
        cds_start = int(fields[6])
        cds_end = int(fields[7])
        if cds_end > cds_start:
            # extract exons in cds range
            exon_number = 0
            for exon_start, exon_end in zip(fields[9].split(','), fields[10].split(',')):
                if exon_start != '':
                    exon_number += 1
                    intersect_range = find_intersect(int(exon_start), int(exon_end), cds_start, cds_end)
                    if intersect_range is not None and intersect_range[1] > intersect_range[0]:
                        item = Interval(start=intersect_range[0], end=intersect_range[1], chrom=chrom)
                        result['cds'].insert(item, other={'name': fields[1], 'strand': fields[3], 'number': exon_number, 'count': int(fields[8])})
                        added += 1
        if i % 10000 == 0:
            write_log(log, 'init_db: {0} lines processed {1} cds intervals last {2}...'.format(i, added, item))
    write_log(log, 'init_db: done with {0} intervals'.format(added))
    return result

def download_db(log):
    '''
        download data from ucsc
    '''
    write_log(log, 'download_db: starting...')
    run("mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e 'select * from refGene' hg19 > gap.db", log)
    write_log(log, 'download_db: done')

def main():
    '''
        parse command line and execute
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Generate gap report')
    parser.add_argument('--min_coverage_ok', required=False, type=int, default=-1, help='maximum value to consider to be low coverage (-1 for all)')
    parser.add_argument('--min_gap_width', required=False, type=int, default=1, help='minimum width of a gap to report')
    parser.add_argument('--coverage', required=True, help='coverage file to examine for gaps')
    parser.add_argument('--db', required=False, help='db to annotate gaps')
    args = parser.parse_args()
    if args.db:
        data_source = init_db(open(args.db, 'r'), sys.stderr)
    else:
        download_db(sys.stderr)
        data_source = init_db(open('gap.db', 'r'), sys.stderr)
    find_gaps(open(args.coverage, 'r'), args.min_gap_width, args.min_coverage_ok, sys.stdout, data_source, sys.stderr)

if __name__ == '__main__':
    main()

