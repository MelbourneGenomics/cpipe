#!/usr/bin/env python
'''
#############################################################################
#
# Melbourne Genomics Coverage Report Script
#
# Copyright Melbourne Genomics Health Alliance members. All rights reserved.
#
# DISTRIBUTION:
#
# This source code should not be distributed to a third party without prior
# approval of the Melbourne Genomics Health Alliance steering committee (via
# Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
#
##############################################################################
#
# Purpose:
#   Generate coverage stats for genes over a capture
# Usage:
#   python calculate_exon_coverage --capture capture_file --exons exon_file
# Outputs:
#   writes a tab separated file of genes and what proportion is covered by the capture
##############################################################################
'''

import collections
import sys

def calculate_coverage(capture, exons, out, log):
    '''
        calculate overlap across genes
    '''
    log.write('reading capture...\n')
    cap = set()
    for line in capture:
        fields = line.strip().split('\t') # chr, start, end
        if len(fields) > 2:
            for x in xrange(int(fields[1]), int(fields[2])):
                cap.add('{0}:{1}'.format(fields[0], x))

    log.write('reading exons...\n')
    found = collections.defaultdict(int)
    total = collections.defaultdict(int)
    for line in exons:
        fields = line.strip().split('\t') # chr, start, end, gene
        if len(fields) > 3:
            gene = fields[3].lower()
            for x in xrange(int(fields[1]), int(fields[2])):
                cand = '{0}:{1}'.format(fields[0], x)
                if cand in cap:
                    found[gene] += 1
                total[gene] += 1

    # write results
    log.write('writing results...\n')
    for gene in sorted(total):
      out.write('{0}\t{1}\n'.format(gene, 100. * found[gene] / total[gene]))

    log.write('done\n')

def main():
    '''
        parse command line and execute
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Generate coverage report')
    parser.add_argument('--capture', required=True, help='capture file')
    parser.add_argument('--exons', required=True, help='exons')
    args = parser.parse_args()
    calculate_coverage(open(args.capture, 'r'), open(args.exons, 'r'), sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()
