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
####################################################################################
#
# Purpose:
#     Given a gene list and reference bed file, generate a bed file with just those genes
# Usage:
#     combine_target_regions.py --genefiles x.gene.txt --bedfiles x.bed --exons exons.bed
# Output files:
#     combination of all regions
####################################################################################
'''

import argparse
import datetime
import os
import sys

def write_log(log, msg):
    '''
        write a date stamped message to log
    '''
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    if log is not None:
        log.write('%s: %s\n' % (now, msg))

def build_genelist(genes, genefiles, beds, log):
    '''
        add genes from genefiles to genes, unless genefile is in beds
    '''
    if genefiles is not None:
        for filename in genefiles:
            if os.path.basename(filename).split('.')[0] in beds:
                write_log(log, 'skipping already specified {0}'.format(filename))
            else:
                if os.path.exists(filename):
                    write_log(log, 'extracting genes from {0}...'.format(filename))
                    added = 0
                    for line in open(filename, 'r'):
                        if line.startswith('#'):
                            continue
                        candidate = line.strip('\n').split('\t')[0]
                        if candidate not in genes:
                            genes.add(candidate.upper().strip())
                            added += 1
                    write_log(log, 'extracting genes from {0}: {1} added'.format(filename, added))
                else:
                    write_log(log, 'WARNING: {0} not found'.format(filename))


def combine_target_regions(genefiles, genefiles_required, bedfiles, exons, target_fh, log):
    '''
        generate a bed file from genefiles and bedfiles
    '''
    # write all bed files unaltered
    beds = set()
    if bedfiles is not None:
        for filename in bedfiles:
            if os.path.exists(filename):
                write_log(log, 'adding bed file {0} to target'.format(filename))
                beds.add(os.path.basename(filename).split('.')[0])
                for line in open(filename, 'r'):
                    target_fh.write(line)
            else:
                write_log(log, 'WARNING: {0} not found'.format(filename))
    # build combined gene list
    genes = set()
    build_genelist(genes, genefiles, beds, log)
    build_genelist(genes, genefiles_required, set(), log)

    # now write filtered exons
    if exons is None:
        if len(genes) > 0:
            write_log(log, 'ERROR: genes specified without exons file')
    else:
        write_log(log, 'writing filtered exons')
        written = 0
        found = set()
        for line in open(exons, 'r'):
            fields = line.strip().split('\t')
            if len(fields) > 3:
                gene = fields[3].strip().upper()
                if gene in genes:
                    found.add(gene)
                    target_fh.write(line)
                    written += 1
        if len(found) == len(genes):
            write_log(log, 'writing filtered exons: {0} genes added in {1} lines'.format(len(found), written))
        else:
            write_log(log, 'WARNING: writing filtered exons: {0} genes extracted from {1} specified. Not found: {2}'.format(len(found), len(genes), ' '.join(sorted(list(genes.difference(found))))))

def main():
    '''
        run from command line
    '''
    parser = argparse.ArgumentParser(description='Generate bed files')
    parser.add_argument('--genefiles', nargs='*', help='list of files containing genes') # input
    parser.add_argument('--genefiles_required', nargs='*', help='list of files containing genes that cannot be skipped') # input
    parser.add_argument('--bedfiles', nargs='*', help='list of bed files to include') # input
    parser.add_argument('--exons', help='exons to convert gene lists to regions')
    parser.add_argument('--skip_same_name', action='store_true', help='skip gene lists with same name as bed file')
    args = parser.parse_args()
    combine_target_regions(args.genefiles, args.genefiles_required, args.bedfiles, args.exons, sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()
