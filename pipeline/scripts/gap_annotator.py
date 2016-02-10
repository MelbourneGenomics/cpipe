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
import sys

# interval matching from https://bitbucket.org/james_taylor/bx-python/raw/ebf9a4b352d3/lib/bx/intervals/operations/quicksect.py
import math
import time
import sys
import random

#class IntervalTree( object ):
#    def __init__( self ):
#        self.chroms = {}
#    def insert( self, interval, linenum=0, other=None ):
#        chrom = interval.chrom
#        start = interval.start
#        end = interval.end
#        if interval.chrom in self.chroms:
#            self.chroms[chrom] = self.chroms[chrom].insert( start, end, linenum, other )
#        else:
#            self.chroms[chrom] = IntervalNode( start, end, linenum, other )
#    def intersect( self, interval, report_func ):
#        chrom = interval.chrom
#        start = interval.start
#        end = interval.end
#        if chrom in self.chroms:
#            self.chroms[chrom].intersect( start, end, report_func )
#    def traverse( self, func ):
#        for item in self.chroms.itervalues():
#            item.traverse( func )
#
#class IntervalNode( object ):
#    def __init__( self, start, end, linenum=0, other=None ):
#        # Python lacks the binomial distribution, so we convert a
#        # uniform into a binomial because it naturally scales with
#        # tree size.  Also, python's uniform is perfect since the
#        # upper limit is not inclusive, which gives us undefined here.
#        self.priority = math.ceil( (-1.0 / math.log(.5)) * math.log( -1.0 / (random.uniform(0,1) - 1)))
#        self.start = start
#        self.end = end
#        self.maxend = self.end
#        self.minend = self.end
#        self.left = None
#        self.right = None
#        self.linenum = linenum
#        self.other = other
#    def insert( self, start, end, linenum=0, other=None ):
#        root = self
#        if start > self.start:
#            # insert to right tree
#            if self.right:
#                self.right = self.right.insert( start, end, linenum, other )
#            else:
#                self.right = IntervalNode(start, end, linenum, other )
#            # rebalance tree
#            if self.priority < self.right.priority:
#                root = self.rotateleft()
#        else:
#            # insert to left tree
#            if self.left:
#                self.left = self.left.insert( start, end, linenum, other )
#            else:
#                self.left = IntervalNode(start, end, linenum, other )
#            # rebalance tree
#            if self.priority < self.left.priority:
#                root = self.rotateright()
#        if root.right and root.left: 
#            root.maxend = max( root.end, root.right.maxend, root.left.maxend )
#            root.minend = min( root.end, root.right.minend, root.left.minend )
#        elif root.right: 
#            root.maxend = max( root.end, root.right.maxend )
#            root.minend = min( root.end, root.right.minend )
#        elif root.left:
#            root.maxend = max( root.end, root.left.maxend )
#            root.minend = min( root.end, root.left.minend )
#        return root
#
#    def rotateright( self ):
#        root = self.left
#        self.left = self.left.right
#        root.right = self
#        if self.right and self.left: 
#            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
#            self.minend = min(self.end, self.right.minend, self.left.minend )
#        elif self.right:
#            self.maxend = max(self.end, self.right.maxend)
#            self.minend = min(self.end, self.right.minend)
#        elif self.left:
#            self.maxend = max(self.end, self.left.maxend)
#            self.minend = min(self.end, self.left.minend )
#        return root
#        
#    def rotateleft( self ):
#        root = self.right
#        self.right = self.right.left
#        root.left = self
#        if self.right and self.left: 
#            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
#            self.minend = min(self.end, self.right.minend, self.left.minend )
#        elif self.right:
#            self.maxend = max(self.end, self.right.maxend)
#            self.minend = min(self.end, self.right.minend)
#        elif self.left:
#            self.maxend = max(self.end, self.left.maxend)
#            self.minend = min(self.end, self.left.minend )
#        return root
#
#    def intersect( self, start, end, report_func ):
#        if start < self.end and end > self.start: report_func( self )
#        if self.left and start < self.left.maxend:
#            self.left.intersect( start, end, report_func )
#        if self.right and end > self.start:
#            self.right.intersect( start, end, report_func )
#
#    def traverse( self, func ):
#        if self.left: self.left.traverse( func )
#        func( self )
#        if self.right: self.right.traverse( func )

def annotate_gap(gap, db):
    '''
        given a gap, annotate it 
    '''
    # find overlapping cds interval, or nearest
    return {}

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

def write_gap(gap, target, db):
    '''
      write out the found gap
    '''
    # note that start and end are inclusive
    annotation = annotate_gap(gap, db)
    target.write('{0},{1},{2},{3},{4},{5},{6},{7}\n'.format(gap['chr'], gap['gene'], gap['start'] + gap['start_offset'] - 1, gap['start'] + gap['start_offset'] + gap['length'] - 1 - 1, min(gap['coverage']), max(gap['coverage']), median(gap['coverage']), gap['length']))

def find_gaps(coverage, min_width, max_coverage, target, db, log):
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
                        write_gap(current, target, db)
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
        write_gap(current, target, db)
    write_log(log, 'finding gaps: {0} lines; {1} gaps: done'.format(idx, gaps))

Interval = collections.namedtuple('Interval', ['start', 'end', 'chrom'])

def init_db(db, log):
    '''
        prepare annotation db
    '''
    write_log(log, 'starting init_db...')
    result = {'cds': IntervalTree()}
    added = 0
    first = True
    for i, line in enumerate(db):
        if first:
            first = False
            continue
        fields = line.strip('\n').split('\t')
        chrom = fields[2]
        cds_start = fields[6].split(',')
        cds_end = fields[7].split(',')
        for j, start in enumerate(cds_start):
            item = Interval(start=int(start), end=int(cds_end[j]), chrom=chrom)
            result['cds'].insert(item)
            added += 1
        if i % 10000 == 0:
            write_log(log, 'init_db: {0} lines processed {1} cds intervals last {2}...'.format(i, added, item))
    write_log(log, 'init_db: done with {0} intervals'.format(added))
    return result

def main():
    '''
        parse command line and execute
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Generate gap report')
    parser.add_argument('--min_coverage_ok', required=False, type=int, default=0, help='maximum value to consider to be low coverage')
    parser.add_argument('--min_gap_width', required=False, type=int, default=1, help='minimum width of a gap to report')
    parser.add_argument('--coverage', required=True, help='coverage file to examine for gaps')
    parser.add_argument('--db', required=False, help='db to annotate gaps')
    args = parser.parse_args()
    if args.db:
        db = init_db(open(args.db, 'r'), sys.stderr)
    else:
        db = {}
    find_gaps(open(args.coverage, 'r'), args.min_gap_width, args.min_coverage_ok, sys.stdout, db, sys.stderr)

if __name__ == '__main__':
    main()
