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
# View the annovarx output more easily
# python ./pipeline/scripts/view_csv.py --line linenum < annovarx.csv
###########################################################################

import argparse
import csv
import sys

def view( fh, num, out ):
  csvfh = csv.reader( fh, delimiter=',', quotechar='"' )
  header = None
  for idx, line in enumerate( csvfh ):
    if not header:
      header = line
    elif num == idx + 1: # 1-based line numbers
      for h, v in zip( header, line ):
        out.write( '{0:>16}: {1}\n'.format( h, v ) )
    else:
      pass
  
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='View CSV file')
  parser.add_argument('--line', required=True, help='line to view')
  args = parser.parse_args()
  view( sys.stdin, int( args.line ), sys.stdout )
