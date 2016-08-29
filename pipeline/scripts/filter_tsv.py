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
#   filter tab separated output from VariantsToTable
# Usage:
#   filter_tsv.py --qual 20 < variants.tsv 1>filtered.tsv
####################################################################################
'''

import argparse
import datetime
import sys

def write_log(log, msg):
    '''
        write a date stamped message to log
    '''
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    if log is not None:
        log.write('%s: %s\n' % (now, msg))

def allow_variant(headers, fields, af_min, qual_min, dp_min, ad_min, stats, proband=''):
    '''
        return true if the variant passes filter cutoffs
    '''
    result = True
    # AF
    try:
        if float(fields[headers.index('AF')]) < af_min:
            result = False
            stats['AF'] += 1
    except ValueError:
        pass

    # QUAL
    try:
        if float(fields[headers.index('QUAL')]) < qual_min:
            result = False
            stats['QUAL'] += 1
    except ValueError:
        pass

    # sample values: DP, AD
    # p is proband; c is count; default DPp to low in case proband is not specified
    sample_results = {'DP': 0, 'DPc': 0, 'DPp': 0}
    passing_samples = set()
    # look at dp
    for header, field in zip(headers, fields):
        if header.endswith('.DP') and field.isdigit():
            # sum of depths
            sample_results['DP'] += int(field)
            sample_results['DPc'] += 1
            # keep track of proband
            if header == '{0}.DP'.format(proband):
                sample_results['DPp'] = int(field)
            # keep track of passing samples
            if int(field) >= dp_min:
                passing_samples.add(header[:-3])

    # look at ad
    ad_ok = False
    for header, field in zip(headers, fields):
        if header.endswith('.AD') and header[:-3] in passing_samples:
            allele_depth = field.split(',')[1]
            if allele_depth.isdigit() and int(allele_depth) >= ad_min:
                ad_ok = True

    if sample_results['DP'] < dp_min * sample_results['DPc'] and sample_results['DPp'] < dp_min: # average dp is too low and proband is too low
        result = False
        stats['DP'] += 1
    if not ad_ok: # ad is too low for all samples with sufficient depth
        result = False
        stats['AD'] += 1
    return result

def filter_tsv(src, target, log, af_min=0.15, qual_min=5, dp_min=5, ad_min=2, proband='', reverse=False):
    '''
        given variants from src, write to target those that pass the filter requirements
    '''
    write_log(log, 'filter_tsv: starting. requirement: af >= {0} qual >= {1} dp >= {2} ad >= {3} proband "{4}"'.format(af_min, qual_min, dp_min, ad_min, proband))
    written = filtered = 0
    stats = {'AF': 0, 'DP': 0, 'QUAL': 0, 'AD': 0}
    header = None
    for line in src:
        if header is None:
            header = line.strip('\n').split('\t')
            target.write(line)
        else:
            fields = line.strip('\n').split('\t')
            # decide whether to keep this line
            if allow_variant(header, fields, af_min, qual_min, dp_min, ad_min, stats, proband=proband):
                if not reverse:
                    target.write(line)
                written += 1
            else:
                if reverse:
                    target.write(line)
                filtered += 1
    total = max(1, written + filtered)
    write_log(log, 'filter_tsv: done. wrote {0} filtered {1} ({2:.1f}%). breakdown: {3}'.format(written, filtered, 100. * filtered / total, ', '.join(['{0}: {1} ({2:.1f}%)'.format(x, stats[x], 100. * stats[x] / total) for x in stats])))

def main():
    '''
        run from command line
    '''
    parser = argparse.ArgumentParser(description='Filter TSVs')
    parser.add_argument('--af', type=float, default=0.15, help='minimum allele frequency')
    parser.add_argument('--qual', type=int, default=5, help='minimum quality')
    parser.add_argument('--dp', type=int, default=5, help='minimum depth')
    parser.add_argument('--ad', type=int, default=2, help='minimum allele depth')
    parser.add_argument('--reverse', action='store_true', default=False, help='write filtered items instead of allowed')
    parser.add_argument('--proband', default='', help='proband sample name')
    args = parser.parse_args()
    filter_tsv(sys.stdin, sys.stdout, sys.stderr, af_min=args.af, qual_min=args.qual, dp_min=args.dp, ad_min=args.ad, proband=args.proband, reverse=args.reverse)

if __name__ == '__main__':
    main()
