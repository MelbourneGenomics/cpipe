#!/usr/bin/env python
####################################################################################
#
# Melbourne Genomics Pipeline Annotation Script
#
# Copyright Melbourne Genomics Health Alliance members. All rights reserved.
#
# DISTRIBUTION:
#
# This source code should not be distributed to a third party without prior
# approval of the Melbourne Genomics Health Alliance steering committee (via
# Natalie Thorne - natalie.thorne@melbournegenomics.org.au).
#
####################################################################################
#
# Purpose:
# * remove spaces intelligently from the priority gene list column
####################################################################################

import re
import sys

GENELIST_COLUMN = 'Prioritised_Genes'

def correct_column( value ):
  '''
    takes spaces out of the gene list while leaving in spaces that separate categories
    @value: uncorrected genelist
    @returns: corrected genelist
  '''
  parts = re.split( '( *[^0-9:]*[0-9]+:)', value.strip() ) # 3:
  corrected = []
  for part in parts:
    if part.endswith(':'):
      corrected.append( re.sub( '^ +', ' ', part ) ) # ensure max 1 space separating gene list
    else:
      corrected.append( part.replace( ' ', '' ) ) # no spaces in actual gene list
  return ''.join( corrected )

def correct_metadata( src, dest ):
  '''
    corrects the genelist column
    @src: file like object containing sample metadata
    @dest: file like object to receive corrected sample metadata
  '''
  first = True
  target_column = None
  for line in src:
    line = line.strip('\n')
    if first:
      fields = line.split( '\t' )
      try:
        target_column = fields.index( GENELIST_COLUMN )
      except ValueError:
        pass
      first = False
      dest.write( '%s\n' % line )
    else:
      if target_column is None:
        dest.write( '%s\n' % line )
      else:
        fields = line.split( '\t' )
        fields[target_column] = correct_column( fields[target_column] )
        dest.write( '%s\n' % ( '\t'.join( fields ) ) )

if __name__ == '__main__':
  correct_metadata( sys.stdin, sys.stdout )

