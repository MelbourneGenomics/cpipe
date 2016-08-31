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
# History:
#   1.0 20-jun-2016 PG initial release
#   1.1 01-jul-2016 PG use bed co-ordinates
#   1.2 31-aug-2016 PG change column ordering, column names, AA calculation
###########################################################################
'''

import collections
import datetime
import gzip
import math
import os
import sys

# interval matching from https://bitbucket.org/james_taylor/bx-python/raw/ebf9a4b352d3/lib/bx/intervals/operations/quicksect.py
import random

__version__ = '1.2'

# 1=include end (compatible style)
# 0=don't include end (bed style)
INCLUDE_END = 0

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

def traversal_handler_builder(gap_start, gap_end):
    '''
        helper to traverse all
    '''
    data = [None, 3e9, None] # current best [interval, distance, direction_to_interval]
    def traversal_handler(interval):
        '''
            closure for traversal
        '''
        if data[0] is None:
            data[0] = interval
            if interval.start >= gap_end:
                data[1] = interval.start - gap_end + 1
                data[2] = 1
            else:
                data[1] = gap_start - interval.end + 1
                data[2] = -1
        elif interval.start >= gap_end and interval.start - gap_end + 1 < data[1]:
            data[0] = interval
            data[1] = interval.start - gap_end + 1
            data[2] = 1
        elif gap_start >= interval.end and gap_start - interval.end + 1 < data[1]:
            data[0] = interval
            data[1] = gap_start - interval.end + 1
            data[2] = -1
    return data, traversal_handler

def annotate_gap_from_ref(gap, data_source, log):
    '''
        given a gap, annotate it
        return annotations, one for each overlap found
    '''
    # find overlapping cds interval, or nearest
    target = []
    if gap['chr'] in data_source['cds'].chroms:
        gap_start_pos = gap['start'] + gap['start_offset'] - 1
        gap_end = gap_start_pos + gap['length']
        data_source['cds'].chroms[gap['chr']].intersect(gap_start_pos, gap_end, lambda x: target.append(x)) # find overlap and add to target
        if len(target) == 0: # no overlap
            # find closest by traversing the entire tree
            result, handler = traversal_handler_builder(gap_start_pos, gap_end)
            data_source['cds'].chroms[gap['chr']].traverse(handler)
            if result[0].other['strand'] == '+':
                exon_rank = result[0].other['number'] # straight from refgene
            else:
                exon_rank = result[0].other['count'] - result[0].other['number'] + 1 # count from the end
            return [{'interval': result[0], 'distance': result[1], 'direction': result[2], 'rank': exon_rank}]
        else: # 1 or more overlaps found
            results = []
            for candidate in target:
                # what is the intersection relative to the coding region?
                abs_intersect = find_intersect(candidate.start, candidate.end, gap_start_pos, gap_end) # absolute position of overlap between exon and gap
                if candidate.other['strand'] == '+': # count from end
                    coding_intersect = [abs_intersect[0] - candidate.start + 1, abs_intersect[1] - candidate.start] # relative position of overlap in the exon. start=1 based, end=inclusive
                else:
                    coding_intersect = [candidate.end - abs_intersect[1] + 1, candidate.end - abs_intersect[0]] # relative position of overlap in the exon. start=1 based, end=inclusive

                #    coding_intersect = [x - target[0].start for x in intersect]
                #else: # count from end
                #    coding_intersect = [target[0].end - x for x in intersect]

                # convert to codon positions i.e. (1, 2, 3, 4, 5) -> (1, 1, 1, 2, 2, ...)
                # TODO strand?
                exon_codon_positions = [int(math.ceil(x / 3.0)) for x in coding_intersect]

                # determine exon rank
                if candidate.other['strand'] == '+':
                    exon_rank = candidate.other['number'] # straight from refgene
                else:
                    exon_rank = candidate.other['count'] - candidate.other['number'] + 1 # count from the end

                # calculate length of all exons preceding this one

                results.append({
                    'interval': candidate,
                    'distance': 0,
                    'coding_intersect': coding_intersect,
                    'exon_codon_positions': exon_codon_positions,
                    'rank': exon_rank})

            #write_log(log, "annotate_gap: coding region intersect: {0} codons: {1}".format(coding_intersect, codon_positions))
            return results
    else:
        return None # unexpected chromosome

def annotate_gap_from_bed(gap, bed, log):
    '''
        find overlaps in specified  bed
    '''
    result = []
    if gap['chr'] in bed.chroms:
        target = []
        gap_start_pos = gap['start'] + gap['start_offset'] - 1
        gap_end = gap_start_pos + gap['length']
        bed.chroms[gap['chr']].intersect(gap_start_pos, gap_end, lambda x: target.append(x)) # find overlap and add to target
        for candidate in target:
            result.append(candidate.other['name'])
    return result

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

def mean(items):
    '''
        geometric average of a set of items
    '''
    if len(items) == 0:
        return 0.

    return 1. * sum(items) / len(items)

def median(items):
    '''
      return the median of a list of items
    '''
    sorted_list = sorted(items)
    if len(items) % 2 == 0:
        high = len(items) / 2
        return (sorted_list[high] + sorted_list[high-1]) / 2.
    else:
        mid = (len(items)- 1) / 2
        return sorted_list[mid]

HEADLINE = ['Chr', 'Start', 'End', 'Gene',
            'Min Cov', 'Max Cov', 'Median Cov', 'Mean Cov', 'Width', 'Tx Name', 'Strand',
            'CDS Distance',
            'CDS Overlap Start', 'CDS Overlap End',
            'Exon Overlap Start', 'Exon Overlap End',
            'AA Overlap Start', 'AA Overlap End',
            'Exon Number', 'Exon Rank']
DEFAULT_NA = 'N/A'

def write_line(target, items):
    '''
        write csv to output
    '''
    target.write('{0}\n'.format(','.join([str(item) for item in items])))

def write_gap(gap, target, data_source, log, beds):
    '''
      write out the found gap
    '''
    # note that we want to write start and end as inclusive

    # annotate from bed files
    additional_headings = []
    additional_data = []
    for bed in sorted(beds.keys()):
        additional_headings.append(bed)
        additional_data.append(';'.join(annotate_gap_from_bed(gap, beds[bed], log))) # join bed

    # annotate from data source
    annotations = annotate_gap_from_ref(gap, data_source, log)
    if annotations is None: # shouldn't happen unless things are really wrong
        write_line(target, [
            gap['chr'],
            gap['start'] + gap['start_offset'] - 1,
            gap['start'] + gap['start_offset'] - 1 + gap['length'] - INCLUDE_END,
            gap['gene'],
            min(gap['coverage']),
            max(gap['coverage']),
            median(gap['coverage']),
            round(mean(gap['coverage']), 1),
            gap['length'],
            DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA] +
                   additional_data)
    else:
        for annotation in annotations: # write a line for each overlap found
            if annotation['distance'] == 0: # gap overlapping coding sequence
                if annotation['interval'].other['strand'] == '+':
                    cds_overlap_start = annotation['interval'].other['previous'] + annotation['coding_intersect'][0]
                    cds_overlap_end = annotation['interval'].other['previous'] + annotation['coding_intersect'][1]
                else:
                    cds_overlap_start = annotation['interval'].other['next'] + annotation['coding_intersect'][0]
                    cds_overlap_end = annotation['interval'].other['next'] + annotation['coding_intersect'][1]
                exon_overlap_start = annotation['coding_intersect'][0]
                exon_overlap_end = annotation['coding_intersect'][1]

                cds_codon_positions = [int(math.ceil(x / 3.0)) for x in (cds_overlap_start, cds_overlap_end)]

                write_line(target, [
                    gap['chr'],
                    gap['start'] + gap['start_offset'] - 1, gap['start'] + gap['start_offset'] - 1 + gap['length'] - INCLUDE_END, gap['gene'],
                    min(gap['coverage']),
                    max(gap['coverage']),
                    median(gap['coverage']),
                    round(mean(gap['coverage']), 1),
                    gap['length'],
                    annotation['interval'].other['name'],
                    annotation['interval'].other['strand'],
                    annotation['distance'],
                    cds_overlap_start, cds_overlap_end,
                    exon_overlap_start, exon_overlap_end,
                    cds_codon_positions[0], cds_codon_positions[1],
                    annotation['interval'].other['number'], annotation['rank']] + additional_data)
            else: # nearest distance
                if annotation['direction'] == 1 and annotation['interval'].other['strand'] == '-' or annotation['direction'] == -1 and annotation['interval'].other['strand'] == '+':
                    distance = -annotation['distance']
                else:
                    distance = annotation['distance']

                write_line(target, [
                    gap['chr'],
                    gap['start'] + gap['start_offset'] - 1,
                    gap['start'] + gap['start_offset'] - 1 + gap['length'] - INCLUDE_END,
                    gap['gene'],
                    min(gap['coverage']),
                    max(gap['coverage']),
                    median(gap['coverage']),
                    round(mean(gap['coverage']), 1),
                    gap['length'],
                    annotation['interval'].other['name'],
                    annotation['interval'].other['strand'],
                    distance,
                    DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA, DEFAULT_NA,
                    annotation['interval'].other['number'],
                    annotation['rank']] + additional_data)

def find_gaps(coverage, min_width, max_coverage, target, data_source, log, beds=None):
    '''
        find gaps and annotate
    '''
    if beds is None:
        beds = {}
    target.write('{0}\n'.format(','.join(HEADLINE + sorted(beds.keys()))))
    current = None
    gaps = 0
    idx = 0
    write_log(log, 'finding gaps...')
    for idx, line in enumerate(coverage):
        fields = line.strip('\n').split('\t') # tab separated: chr, start, end, gene, offset, coverage
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
                        write_gap(current, target, data_source, log, beds)
                        gaps += 1
                    if max_coverage == -1 or coverage <= max_coverage: # start a new gap?
                        current = {
                            'start': int(fields[1]),
                            'start_offset': int(fields[4]),
                            'length': 1,
                            'chr': fields[0],
                            'gene': fields[3],
                            'coverage': [coverage]}
                    else:
                        current = None
        else:
            write_log(log, 'skipped line {0}'.format(idx))

        if idx % 100000 == 0:
            write_log(log, 'finding gaps: {0} lines processed; {1} gaps found...'.format(idx, gaps))

    # still a gap in progress?
    if current is not None and current['length'] >= min_width:
        write_gap(current, target, data_source, log, beds)
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

def init_db(target, log, exclusions=None):
    '''
        prepare annotation db
        @target: is refgene, and expects the following fields:
         - name(1), chrom(2), strand(3), cds_start(6), cds_end(7), count(8), exon_start(9), exon_end(10)
    '''
    write_log(log, 'starting init_db...')
    if exclusions is None:
        exclusions = set()
    result = {'cds': IntervalTree()}
    added = 0
    first = True
    position = 0
    for position, line in enumerate(target):
        if first:
            first = False
            continue
        fields = line.strip('\n').split('\t')
        name = fields[1]
        if name.strip().upper() in exclusions:
            continue

        chrom = fields[2]
        cds_start = int(fields[6])
        cds_end = int(fields[7])
        if cds_end > cds_start:
            # extract exons in cds range
            exon_number = 0
            exon_sum = 0
            to_add = []
            sizes = []
            strand = fields[3]
            for exon_start, exon_end in zip(fields[9].split(','), fields[10].split(',')):
                if exon_start != '':
                    exon_number += 1
                    intersect_range = find_intersect(int(exon_start), int(exon_end), cds_start, cds_end)
                    if intersect_range is not None and intersect_range[1] > intersect_range[0]:
                        item = Interval(start=intersect_range[0], end=intersect_range[1], chrom=chrom)
                        other = {
                            'name': name,
                            'strand': strand,
                            'number': exon_number,
                            'count': int(fields[8]),
                            'previous': exon_sum}
                        to_add.append({'item': item, 'other': other}) # previous contains sum of exons before this one
                        sizes.append(intersect_range[1] - intersect_range[0])
                        exon_sum += sizes[-1]
                        added += 1
            # calculate upcoming exon sum
            to_add[-1]['other']['next'] = 0
            for i, exon in enumerate(to_add):
                exon['other']['next'] = sum(sizes[i+1:]) # next contains sum of exons after this one
                result['cds'].insert(exon['item'], other=exon['other'])
        if position % 10000 == 0:
            write_log(log, 'init_db: {0} lines processed {1} cds intervals last {2}...'.format(position, added, item))
    write_log(log, 'init_db: done with {0} intervals'.format(added))
    return result

def init_bed(target, target_name, log):
    '''
        prepare interval file for bed
    '''
    write_log(log, 'processing {0}...'.format(target_name))
    result = IntervalTree()
    i = 0
    for i, line in enumerate(target):
        fields = line.strip('\n').split('\t')
        if len(fields) < 4:
            write_log(log, 'WARNING: {0} contains too few columns on line {1}'.format(target_name, i))
        result.insert(Interval(start=int(fields[1]), end=int(fields[2]), chrom=fields[0]), other={'name': fields[3]})
    write_log(log, 'processing {0}: done processing {1} lines'.format(target_name, i))
    return result

def download_db(log):
    '''
        download data from ucsc
    '''
    write_log(log, 'download_db: starting...')
    run("mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e 'select * from refGene' hg19 > gap.db", log)
    write_log(log, 'download_db: done')

def build_exclusion_list(incoming_fh):
    '''
        build a list of transcripts to exclude
    '''
    result = set()
    for line in incoming_fh:
        result.add(line.strip().upper())
    return result

def main():
    '''
        parse command line and execute
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Generate gap report')
    parser.add_argument('--max_low_coverage', required=False, type=int, default=-1, help='maximum value to consider to be low coverage (-1 for all)')
    parser.add_argument('--min_gap_width', required=False, type=int, default=1, help='minimum width of a gap to report')
    parser.add_argument('--coverage', required=True, help='coverage file to examine for gaps')
    parser.add_argument('--db', required=False, help='db to annotate gaps')
    parser.add_argument('--exclude', required=False, help='file containing transcripts to exclude from gap report')
    parser.add_argument('--beds', required=False, nargs='*', help='additional bed files to find overlaps in')
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {0}'.format(__version__))
    args = parser.parse_args()

    # exclusions
    if args.exclude:
        transcript_exclusions = build_exclusion_list(open(args.exclude, 'r'))
    else:
        transcript_exclusions = set()

    # refgene
    if args.db:
        data_source = init_db(open(args.db, 'r'), sys.stderr, exclusions=transcript_exclusions)
    else:
        download_db(sys.stderr)
        data_source = init_db(open('gap.db', 'r'), sys.stderr, exclusions=transcript_exclusions)

    # bed files
    beds = {}
    if args.beds:
        for bed in args.beds:
            target_name = os.path.basename(bed).split('.')[0]
            beds[target_name] = init_bed(open(bed, 'r'), target_name, log=sys.stderr)

    # coverage data
    if args.coverage.endswith('.gz'):
        coverage_fh = gzip.open(args.coverage, 'r')
    else:
        coverage_fh = open(args.coverage, 'r')

    find_gaps(coverage_fh, args.min_gap_width, args.max_low_coverage, sys.stdout, data_source, sys.stderr, beds)

if __name__ == '__main__':
    main()

