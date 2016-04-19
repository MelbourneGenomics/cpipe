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
#   Determine where variants were filtered
#   Generates a tsv of the form: sample\t chrom\t pos\t seen
#     where seen is exome, cohort, or FINAL
# Usage:
#   See main() for arguments
#
###########################################################################
'''

import collections
import glob
import os.path
import sys

def generate_report(source_dir, log):
    '''
        generate data showing which variants are in which vcf files
    '''
    log.write('filter_report: reading files from {0}\n'.format(source_dir))
    data = collections.defaultdict(dict)
    processed = 0
    for filename in glob.glob(os.path.join(source_dir, '*.vcf')):
        log.write('processing {0}...\n'.format(filename))
        base = os.path.basename(filename)
        sample = base.split('.')[0]
        if sample not in data:
            data[sample] = collections.defaultdict(set)
        stage = '.'.join(base.split('.')[1:-1])
        for line in open(base, 'r'):
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) > 2:
                key = '{0}\t{1}'.format(fields[0], fields[1])
                data[sample][key].add(stage)
                processed += 1
                if processed % 100000 == 0:
                    log.write('processed {0} records\n'.format(processed))
    log.write('filter_report: done processing {0} records\n'.format(processed))
    return data

def to_seen(seen):
    '''
        string representation of set of stages
        assume stages are:
            merge.dedup.realign.recal.g
            merge.dedup.realign.recal.g.filter_variants.merge_variants_gvcf
            merge.dedup.realign.recal.genotype.raw
            genotype.norm
            genotype.soi
            genotype.soi.vep
            genotype.soi.vep.post_filter
    '''
    if 'genotype.soi.vep.post_filter' in seen:
        return 'FINAL'
    elif 'merge.dedup.realign.recal.genotype.raw' in seen:
        return 'genotype filter'
    elif 'merge.dedup.realign.recal.g.filter_variants.merge_variants_gvcf' in seen:
        return 'cohort'
    else:
        return 'exome'
    
    
def write_report(data, target, log):
    '''
        write tsv to target
    '''
    log.write('write_report: starting...\n')
    target.write('{0}\t{1}\t{2}\t{3}\n'.format('Sample', 'Chr', 'Pos', 'Seen'))
    for sample in sorted(data.keys()):
        log.write('write_report: sample {0}...\n'.format(sample))
        for variant in sorted(data[sample].keys()):
            #sys.stderr.write('sample {0} variant {1}\n'.format(sample, variant))
            target.write('{0}\t{1}\t{2}\t{3}\n'.format(sample, variant.split('\t')[0], variant.split('\t')[1], to_seen(data[sample][variant])))
            
    log.write('write_report: done\n')

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Variant filter report')
    parser.add_argument('--source_dir', default='.', help='Get VCFs from here')
    args = parser.parse_args()
    report = generate_report(source_dir=args.source_dir, log=sys.stderr)
    write_report(data=report, target=sys.stdout, log=sys.stderr)

if __name__ == '__main__':
    main()
