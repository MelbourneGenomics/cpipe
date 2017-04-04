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
####################################################################################
#
# Purpose:
# * Find genes missing from a bed file
#
# Usage:
#   find_missing_genes bed_file < genes
####################################################################################

import sys

def main():
  # read bed genes
  ref = set()
  for line in open( sys.argv[1], 'r' ):
    fields = line.strip().split( '\t' )
    if len(fields) > 3:
      ref.add( fields[3].upper() )

  # read genes
  missing = set()
  for line in sys.stdin:
    if line.startswith( '#' ):
      continue
    fields = line.strip().split('\t')
    candidate = fields[0].upper()
    if candidate not in ref and candidate != '1' and candidate != 'HGNC_SYMBOL':
      missing.add( candidate )

  print(('\n'.join( sorted( list( missing ) ) )))

if __name__ == '__main__':
  main()
