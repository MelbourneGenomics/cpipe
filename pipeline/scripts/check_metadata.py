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
# This script takes a sample metadata file as input (stdin) and
# prints out the field values in a human readable form
###########################################################################
'''

import sys

is_numeric = set( [ 'dna concentration', 'dna quantity', 'dna quality', 'mean coverage' ] )
is_enumeration = { 'sex': set( [ 'Male', 'Female', 'Unknown', 'other' ] ), 'sample type': set( [ 'Normal', 'Tumour' ] ), 'consanguinity': set( [ 'No', 'Yes', 'Suspected', 'Unknown' ] ), 'ethnicity': set( [ 'Unknown', 'European', 'African', 'Asian' ] ) }
is_date = set( [ 'dna_date', 'capture_date', 'sequencing_date' ] )

def is_valid_numeric( field ):
  if len(field) > 0:
    try:
      float(field)
      return True
    except ValueError:
      return False
  return True

def is_valid_enumeration( field, allowed ):
  return len(field) == 0 or field in allowed

def is_valid_date( field ):
  return len(field) == 0 or len(field) == 8 and field.isdigit()

def validate( fh, out, err ):
  headers = fh.readline().strip().split('\t')
  idx = -1
  warnings = []
  for idx, line in enumerate(fh):
    out.write( "===== Sample {0} =====\n".format( idx ) )
    fields = line.strip('\n').split('\t')
    for jdx, field in enumerate(fields):
      out.write( "{0:>24}: {1}\n".format( headers[jdx], field ) )
      if field.startswith( ' ' ):
        warnings.append( 'Sample {0} field {1} (column {2}) contains leading whitespace'.format( idx, headers[jdx], jdx ) )
      if field.endswith( ' ' ):
        warnings.append( 'Sample {0} field {1} (column {2}) contains trailing whitespace'.format( idx, headers[jdx], jdx ) )
      if headers[jdx].lower() in is_numeric and not is_valid_numeric( field ):
          warnings.append( 'Sample {0} field {1} (column {2}) cannot be "{3}": must be empty or a number'.format( idx, headers[jdx], jdx, field ) )
      if headers[jdx].lower() in is_enumeration and not is_valid_enumeration( field, is_enumeration[headers[jdx].lower()] ):
          warnings.append( 'Sample {0} field {1} (column {2}) cannot be "{3}": must be empty or one of: {4}'.format( idx, headers[jdx], jdx, field, ', '.join( is_enumeration[headers[jdx].lower()] ) ) )
      if headers[jdx].lower() in is_date and not is_valid_date( field ):
          warnings.append( 'Sample {0} field {1} (column {2}) cannot be "{3}": must be empty or date (yyyymmdd)'.format( idx, headers[jdx], jdx, field ) )

  if idx == -1:
    err.write( "ERROR: file only contains one line. Are you using Windows style line feeds?\n" )
  for warning in warnings:
    err.write( "WARNING: {0}\n".format( warning ) )
  if len(warnings) == 0:
    err.write( "No warnings\n" )

if __name__ == '__main__':
  validate( sys.stdin, sys.stdout, sys.stderr )
