#!/usr/bin/env python
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

import os
import csv
from subprocess import call
from argparse import (ArgumentParser, FileType, ArgumentDefaultsHelpFormatter)

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Produce a bam file of reads overlaping a variant. By default, reads overlaping the region 100bp upstream and downstream of the variant are included.',
    formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--bam', type=str, required=True,
        help='Input bam file.')
    parser.add_argument(
        '--tsv', type=str, required=True,
        help='Variants in tsv format. Must have at least the columns CHROM, POS, REF, ALT, HGVSc. Sample assumed to be the start of this file name (before the ".").')
    parser.add_argument(
        '--outdir', type=str, default='.',
        help='Output directory for bams. Directory must already exist.')
    parser.add_argument(
        '--log', type=str,
        help='Optionally output the number of bams created to a log file with this name.')
    parser.add_argument(
        '--upstream', type=int, default=100,
        help="Include reads this many bp upstream (5') of the variant.")
    parser.add_argument(
        '--downstream', type=int, default=100,
        help="Include reads this many bp downstream (3') of the variant.")
    parser.add_argument(
        '--samtoolsdir', type=str,
        help="Directory in which samtools is installed. Default: samtools assumed to be in PATH.")

    return parser.parse_args() 


def main():
    # Parse command line arguments
    args = parse_args()
    inbam = args.bam
    variantfile = args.tsv
    outdir = args.outdir
    upstream = args.upstream # Default 100 
    downstream = args.downstream # Default 100
    if args.samtoolsdir:
        samtools_exec = args.samtoolsdir + '/samtools'
    else:
        samtools_exec = 'samtools'

    # Assume sample name is the first part of the filename before the "."
    #sample = variantfile.split('/')[-1].split('.')[0]
    variant_filename = variantfile.split('/')[-1] # vlsci_1.0.2_000000017_012345678.lovd.tsv
    if '_' in variant_filename:
      # take the last part after the _ and the first part before the .
      #sample = variant_filename.split('_')[-1].split('.')[0]
      sample = '.'.join( variant_filename.split('.')[:-2] ) # take everything before the last 2 dots
    else:
      # Assume sample name is the first part of the filename before the "."
      sample = variant_filename.split('.')[0]

    with open(variantfile) as varianttsv:
        var_count = 0 
        for line in csv.DictReader(varianttsv, delimiter='\t'):
            # annovar NM = line['AAChange'].split(':')[0]
            nm = line['HGVSc'].split(':')[0].replace('.', '_')

            chromosome = line['CHROM']
            start = int(line['POS']) # annovar int(line['Start'])
            variation_length = len(line['ALT']) - len(line['REF'])
            if variation_length > 0:
                end = start + variation_length # annovar int(line['End'])
            else:
                end = start

            # annovar outbam = '{0}-{1}-{2}-{3}-{4}-IGV.bam'.format(sample, NM, chr, start, end)
            outbam = '{0}-{1}-{2}-{3}-{4}-IGV.bam'.format(sample, nm, chromosome, start, end)
            region = '{0}:{1}-{2}'.format(chromosome, start - upstream, end + downstream)

            # create new bam containing reads in the given region
            call([samtools_exec, 'view', '-b', '-o', outdir+'/'+outbam, inbam, region])
            # index the bam
            call([samtools_exec, 'index', outdir+'/'+outbam])
            
            var_count += 1

    # If required, produce a log file
    # This is useful to trick bpipe into tracking the output of this script, 
    # since an unknown number of bam files are produced (sometimes zero)
    if args.log:
        with open(args.log, 'w') as logfile:
            logfile.write('{0} bams for sample {1} written to {2}'.format(var_count, sample, outdir))

if __name__ == '__main__':
    main()
