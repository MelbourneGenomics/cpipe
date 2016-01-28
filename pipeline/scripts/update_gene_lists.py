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
###################################################################################
#
# Purpose:
#   Given a gene list and reference bed file, generate a bed file with just those genes
#   The source file is source_dir/cohort.genes.txt
#   The target file is target_dir/cohort/cohort.genes.txt
# Usage:
#   update_gene_lists --source dir --target dir
####################################################################################
'''

import argparse
import datetime
import glob
import os
import sys

CATEGORY = '1'

def write_log(log, msg):
    '''
        write date stamped msg to log
    '''
    log.write('%s: %s\n' % (datetime.datetime.now().strftime('%y%m%d-%H%M%S'), msg))

def update_gene_lists(source_dir, target_dir, log):
    '''
        adds genes from files of the form source_dir/*.add.genes.txt to gene lists in target_dir/cohort/cohort.genes.txt
    '''
    for filename in glob.glob(os.path.join(source_dir, '*.add.genes.txt')):
        cohort = os.path.basename(filename).split('.')[0]
        # find corresponding flagship
        target = os.path.join(target_dir, cohort, '%s.genes.txt' % cohort) # target/cohort/cohort.genes.txt
        if os.path.isfile(target):
            # read existing genes and categories
            genes = {}
            for line in open(target, 'r'):
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                genes[fields[0].upper()] = fields[1]
            # read new genes
            added = set()
            candidates = set()
            for line in open(filename, 'r'):
                if line.startswith('#'):
                    continue
                gene = line.strip().upper()
                candidates.add(gene)
                if gene not in genes:
                    added.add(gene)
                    genes[gene] = CATEGORY

            # write out additional
            if len(added) > 0:
                with open(target, 'w') as fh:
                    fh.write('#version %s\n' % datetime.datetime.now().strftime('%y%m%d'))
                    fh.write('#notes %i gene(s) added: %s\n' % (len(added), ','.join(sorted(list(added)))))
                    for gene in sorted(genes.keys()):
                        fh.write('%s\t%s\n' % (gene, genes[gene]))
                write_log(log, '%s: %i gene(s) added from %i candidate(s): %s' % (cohort, len(added), len(candidates), ','.join(sorted(list(added)))))
            else:
                write_log(log, '%s: no changes from %i candidate(s)' % (cohort, len(candidates)))
        else:
            write_log(log, 'ERROR: target gene list %s does not exist' % target)

def main():
    '''
        update gene lists from command line options
    '''
    parser = argparse.ArgumentParser(description='Generate bed files')
    parser.add_argument('--source', required=True, help='source of extra genes') # input
    parser.add_argument('--target', required=True, help='target containing gene files to update') # input
    parser.add_argument('--log', required=False, help='write changes to this file') # input
    args = parser.parse_args()
    if args.log:
        log = open(args.log, 'a+')
    else:
        log = sys.stderr

    update_gene_lists(args.source, args.target, log)

if __name__ == '__main__':
    main()
