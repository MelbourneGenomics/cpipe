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
# Extracts information about variants found for a particular gene
# e.g.
# python ./pipeline/scripts/examine_variant.py --gene ACTN2 --batch na18507.150907 --sample NA18507
###########################################################################

import argparse
import csv
import glob
import sys

def report_gene( name, gene, fn, out, separator, display, use_header=False ):
  out.write( '===== {0} =====\n'.format( name ) )
  found = -1
  for file in glob.glob( fn ):
    first = True
    with open( file, 'r' ) as fh:
      out.write( '{0}\n'.format( file ) )
      csvfh = csv.reader( fh, delimiter=separator, quotechar='"' )
      for line in csvfh:
        if first and use_header:
          header = line
        if first or any( [ gene in cell for cell in line ] ):
          first = False
          found += 1
          filtered = []
          for idx, field in enumerate(line):
            if use_header and idx < len(header):
              key = header[idx]
            else:
              key = idx
            if key in display:
              filtered.append( field )
          out.write( '{0}: {1}\n'.format( gene, '\t'.join( filtered ) ) )
  out.write( '--- {0} instances found ---\n'.format( found ) )

def examine( gene, batch, sample, out ):
  # stage 1 is VEP: *NA18507.merge.dedup.realign.recal.filter_variants.merge_variants.vep.sort.vcf  
  report_gene( 'VEP', gene, './batches/{0}/analysis/variants/*{1}.merge.dedup.realign.recal.filter*.vep.sort.vcf'.format( batch, sample ), sys.stdout, separator='\t', display=set( [0, 1, 3, 4] ) )
  report_gene( 'Annovar', gene, './batches/{0}/analysis/variants/*{1}.merge.dedup.realign.recal.filter*.vep.sort.hg19_multianno.csv'.format( batch, sample ), sys.stdout, separator=',', display=set( ['Chr', 'Start', 'Ref', 'Alt', 'Func', 'Priority_Index'] ), use_header=True )
  report_gene( 'Significance', gene, './batches/{0}/analysis/variants/*{1}.merge.dedup.realign.recal.filter*.vep.sort.hg19_multianno.con.sig.csv'.format( batch, sample ), sys.stdout, separator=',', display=set( ['Chr', 'Start', 'Ref', 'Alt', 'Func', 'Priority_Index'] ), use_header=True )
  report_gene( 'Final', gene, './batches/{0}/analysis/results/*{1}.annovarx.csv'.format( batch, sample ), sys.stdout, separator=',', display=set( ['Chr', 'Start', 'Ref', 'Alt', 'Priority_Index', '#Obs'] ), use_header=True )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Examine variant filtering')
  parser.add_argument('--gene', required=True, help='name of gene to examine')
  parser.add_argument('--sample', required=True, help='sample to examine')
  parser.add_argument('--batch', required=True, help='batch to examine')
  args = parser.parse_args() 
  examine( args.gene, args.batch, args.sample, sys.stdout )
