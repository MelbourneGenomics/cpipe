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
###################################################################################
#
# Purpose:
#   Converts Annovar's hg19_refGene.txt to a bed file
#   Use in conjunction with refgene_to_bed.sh
#   Not used as part of the pipeline run, this is an administrative tool.
#
####################################################################################
import datetime
import sys

if len(sys.argv) > 1 and sys.argv[1] == 'post': # split out genes
  sys.stdout.write( '#version %s\n' % datetime.datetime.now().strftime("%Y%m%d") )
  for line in sys.stdin:
    if line.startswith( '#' ):
      sys.stdout.write( line )
    else:
      fields = line.strip().split('\t')
      genes = set( fields[3].split(';') )
      for gene in genes:
        sys.stdout.write( '%s\t%s\t%s\t%s\n' % ( fields[0], fields[1], fields[2], gene ) )
else: # convert from refGene
  sys.stdout.write( '#version %s\n' % datetime.datetime.now().strftime("%Y%m%d") )
  for line in sys.stdin:
    fields = line.strip().split( '\t' )
    exons = ( fields[9].split(','), fields[10].split(',') )
    for exon_start, exon_end in zip( *exons ):
      if exon_start != '' and exon_end != '':
        sys.stdout.write( '%s\t%s\t%s\t%s\n' % ( fields[2], exon_start, exon_end, fields[12] ) )
