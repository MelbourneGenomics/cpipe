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
# Compare two analyses.
# This looks at the annovarx.csv output and shows differences between the two
# python ~/compare_analyses.py --dir1 ./prod/batches/b1 --dir2 ./stage/batches/b2 --sample1 123456789 --sample2 9877654321
#
# In the future we may look at previous stages to determine why a variant was filtered
# In the meantime, look at evaluate_variant.py
###########################################################################

import argparse
import csv
import glob
import sys

def find_variants( fh ):
  csvfh = csv.reader( fh, delimiter=',', quotechar='"' )
  indexes = None
  data_indexes = None
  result = set()
  extra = {}
  for line in csvfh:
    if not indexes: # first line
      indexes = [ line.index(x) for x in ('Gene','Chr','Start') ]
      data_indexes = [ line.index(x) for x in ('Gene','Chr','Start', 'Func') ]
    else:
      key = [ line[i] for i in indexes ]
      data = [ line[i] for i in data_indexes ]
      # hack to deal with ;
      if ';' in key[0]:
        for i, gene in enumerate( key[0].split(';') ):
          skey = '\t'.join( [ gene, key[1], key[2] ] )
          data_fixed = '\t'.join( [ gene, data[1], data[2], data[3].split(';')[i] ] )
          result.add( skey )
          extra[ skey ] = data_fixed
      else:
        skey = '\t'.join( key )
        result.add( skey )
        extra[ skey ] = '\t'.join( data )
  return result, extra
  
def compare( d1, d2, s1, s2, out, common=False ):
  # compare the annovars
  a1fn = glob.glob( '{0}/analysis/results/*{1}.annovarx.csv'.format( d1, s1 ) )[0]
  a2fn = glob.glob( '{0}/analysis/results/*{1}.annovarx.csv'.format( d2, s2 ) )[0]
  a1, a1extra = find_variants( open( a1fn, 'r' ) )
  a2, a2extra = find_variants( open( a2fn, 'r' ) )
  out.write( '{0} total variants in {1} {2}\n'.format( len(a1), d1, s1 ) )
  out.write( '{0} total variants in {1} {2}\n'.format( len(a2), d2, s2 ) )
  # common
  if common:
    both = a1.intersection(a2)
    out.write( '----- {0} variants in common -----\n'.format( len(both) ) )
    for x in sorted( list( both ) ):
      out.write( '{0}\n'.format( a1extra[x] ) )
  # only s1
  s1only = a1.difference( a2 )
  out.write( '----- {0} variants only in {1} {2} -----\n'.format( len(s1only), d1, s1 ) )
  for x in sorted( list( s1only ) ):
    out.write( '{0}\n'.format( a1extra[x] ) )
  # only s2
  s2only = a2.difference( a1 )
  out.write( '----- {0} variants only in {1} {2} -----\n'.format( len(s2only), d2, s2 ) )
  for x in sorted( list( s2only ) ):
    out.write( '{0}\n'.format( a2extra[x] ) )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Compare two analyses')
  parser.add_argument('--dir1', required=True, help='batch 1 directory')
  parser.add_argument('--dir2', required=True, help='batch 2 directory')
  parser.add_argument('--sample1', required=True, help='sample 1 name')
  parser.add_argument('--sample2', required=True, help='sample 2 name')
  parser.add_argument('--common', action='store_true', required=False, default=False, help='show variants in common' )
  args = parser.parse_args()
  compare( args.dir1, args.dir2, args.sample1, args.sample2, sys.stdout, common=args.common )
