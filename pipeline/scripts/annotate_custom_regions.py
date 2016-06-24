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
#   Adds a field name to a TSV indicating if any variant overlaps with specified regions
# Usage:
#   annotate_custom_regions --bed custom.bed < before.tsv > after.tsv
#
# arguments:
# --bed: bed file to check
#
##############################################################################
'''

import argparse
import collections
import math
import sys

FIELD_NAME = 'CPIPE_BED'
NOT_FOUND = ''
FOUND = '1'
CHROM = 'CHROM'
POS = 'POS'
SEPARATOR = ';' # if found more than once

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

Interval = collections.namedtuple('Interval', ['start', 'end', 'chrom'])

def find_region(chrom, position, regions):
    '''
        returns the name of the overlap else None if no overlap found
    '''
    result = None
    if chrom in regions.chroms:
        target = []
        regions.chroms[chrom].intersect(position, position + 1, lambda x: target.append(x)) # find overlap and add to target
        if len(target) > 0:
            return SEPARATOR.join(sorted(set([x.other['name'] for x in target])))
    return result

def annotate(source, target, regions, log):
    '''
        add annotation field to target by reading source and considering regions
    '''
    log.write('annotate: starting...\n')
    first = True
    found = 0
    for line in source:
        if first:
            first = False
            header = line.strip('\n').split('\t')
            if CHROM not in header or POS not in header:
                log.write('ERROR: fields {0} and {1} not found in {2}\n'.format(CHROM, POS, header))
                sys.exit(1)
            header.append(FIELD_NAME)
            target.write('{0}\n'.format('\t'.join(header)))
            continue
        fields = line.strip('\n').split('\t')
        # extract chromosome and position
        chrom = fields[header.index(CHROM)]
        position = int(fields[header.index(POS)])
        # is it in a region?
        name = find_region(chrom, position, regions)
        if name is None:
            fields.append(NOT_FOUND)
        else:
            fields.append(name)
            found += 1
        target.write('{0}\n'.format('\t'.join(fields)))
    log.write('annotate: found {0} overlapping positions\n'.format(found))

def build_regions(bed_fh, log):
    '''
        build the interval tree from the bed file
    '''
    log.write('processing bed file...\n')
    result = IntervalTree()
    lines = 0
    for lines, line in enumerate(bed_fh):
        fields = line.strip('\n').split('\t')
        if len(fields) < 3:
            log.write('WARNING: bed file contains too few columns on line {1}\n'.format(lines))
        if len(fields) < 4:
            result.insert(Interval(start=int(fields[1]), end=int(fields[2]), chrom=fields[0]), other={'name': FOUND})
        else:
            result.insert(Interval(start=int(fields[1]), end=int(fields[2]), chrom=fields[0]), other={'name': fields[3]})
    log.write('processing bed file: done processing {0} lines...\n'.format(lines + 1))
    return result

def main(source, target, log):
    '''
        run from command line
    '''
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--bed', required=True, help='bed file containing regions to search')
    args = parser.parse_args()
    annotate(source, target, build_regions(open(args.bed, 'r'), log), log)

if __name__ == '__main__':
    main(sys.stdin, sys.stdout, sys.stderr)

