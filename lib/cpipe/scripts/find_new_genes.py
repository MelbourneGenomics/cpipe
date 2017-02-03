#!/usr/bin/env python
"""
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
#     find_new_genes --reference reference_bed --exclude exclude --target target_dir < sample_file
# Output files:
#     - cohort.add.genes.txt
#     - cohort.addonce.genes.txt
#     - cohort.notfound.genes.txt
#     - cohort.excluded.genes.txt
####################################################################################
"""

import argparse
import os
import re
import sys

def generate_new_genes(sample_lines, log, reference_genes, excluded_genes, reference_source, excluded_source):
    """
        given samples, make files of the form CS.extra.genes.txt and CS.extra.excluded.genes.txt
    """
    # get list of reference genes
    reference = set()
    for line in reference_genes:
        fields = line.strip().split('\t')
        if len(fields) > 3:
            reference.add(fields[3].upper())
    # get list of excluded genes
    excluded = set()
    for line in excluded_genes:
        if line.startswith('#'):
            continue
        excluded.add(line.strip().split('\t')[0].upper())
    log.write('{0} available reference genes found in {1}, {2} excluded genes found in {3}: {4}\n'.format(len(reference), reference_source, len(excluded), excluded_source, ' '.join(sorted(list(excluded)))))

    # parse sample
    headers = {}
    for idx, title in enumerate(sample_lines[0].split('\t')):
        headers[title] = idx

    if 'Cohort' not in headers:
        log.write('ERROR: Cohort not in sample header\n')
        return 1
    if 'Prioritised_Genes' not in headers:
        log.write('ERROR: Prioritised_Genes not in sample header\n')
        return 1

    result = {}
    for line in sample_lines[1:]: # each sample
        fields = line.split('\t')
        cohort = fields[headers['Cohort']]
        genes = fields[headers['Prioritised_Genes']]
        sample_id = fields[headers['Sample_ID']]
        addonce_target = 'addonce.{0}'.format(sample_id)
        if cohort != '':
            if cohort not in result:
                result[cohort] = {'add': set(), 'notfound': set()}
            result[cohort][addonce_target] = set() # need an addonce for every sample (even if it's empty)
            candidates = re.split('[:,"]+', genes.upper())
            for candidate in candidates:
                if candidate in excluded:
                    result[cohort][addonce_target].add(candidate)
                elif candidate in reference:
                    result[cohort]['add'].add(candidate)
                else:
                    if candidate != '' and candidate != '4':
                        result[cohort]['notfound'].add(candidate)
    return result

def write_genes(additions, target, log, dummy=False):
    """
        writes e.g. CS.add.genes.txt, CS.excluded.genes.txt, CS.notfound.genes.txt
    """
    for cohort in additions:
        for category in additions[cohort]:
            if dummy:
                for gene in additions[cohort][category]:
                    log.write('%s.%s.%s\n' % (cohort, category, gene))
            else:
                with open(os.path.join(target, '%s.%s.genes.txt' % (cohort, category)), 'w') as fh_out:
                    for gene in additions[cohort][category]:
                        fh_out.write('%s\n' % gene)

def main():
    """
        run from command line
    """
    parser = argparse.ArgumentParser(description='Generate bed files')
    parser.add_argument('--reference', required=True, help='reference bed file') # input
    parser.add_argument('--exclude', required=True, help='file containing genes to exclude') # input
    parser.add_argument('--target', required=True, help='target directory')
    args = parser.parse_args()
    samples = sys.stdin.readlines()
    additions = generate_new_genes(samples, sys.stderr, open(args.reference, 'r'), open(args.exclude, 'r'), os.path.basename(args.reference), os.path.basename(args.exclude))
    write_genes(additions, args.target, sys.stderr, dummy=False)

if __name__ == '__main__':
    main()
