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
# convert a bed file into a list of genes
############################################################################

import sys

genes = set()
for line in sys.stdin:
  if line.startswith( '#' ):
    continue
  fields = line.strip().split()
  if len(fields) > 3:
    gene = fields[3].strip().upper()
    genes.add( gene )

for gene in sorted( list( genes ) ):
  sys.stdout.write( '{0}\t{1}\n'.format( gene, 1 ) )
