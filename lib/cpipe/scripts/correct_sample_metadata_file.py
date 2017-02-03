#!/usr/bin/env python
"""
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
# * remove spaces intelligently from the priority gene list column
####################################################################################
"""

import re
import sys

GENELIST_COLUMN = 'Prioritised_Genes'

def correct_column(value):
    """
        takes spaces out of the gene list while leaving in spaces that separate categories
        @value: uncorrected genelist
        @returns: corrected genelist
    """
    parts = re.split('( *[^0-9:]*[0-9]+:)', value.strip()) # 3:
    corrected = []
    for part in parts:
        if part.endswith(':'):
            corrected.append(re.sub('^ +', ' ', part)) # ensure max 1 space separating gene list
        else:
            corrected.append(part.replace(' ', '')) # no spaces in actual gene list
    return ''.join(corrected)

def correct_metadata(src, dest):
    """
        corrects the genelist column
        @src: file like object containing sample metadata
        @dest: file like object to receive corrected sample metadata
    """
    first = True
    target_column = None
    for line in src:
        line = line.strip('\n')
        if first:
            fields = line.split('\t')
            try:
                target_column = fields.index(GENELIST_COLUMN)
            except ValueError:
                pass
            first = False
            dest.write('%s\n' % line)
        else:
            if target_column is None:
                dest.write('%s\n' % line)
            else:
                fields = line.split('\t')
                fields[target_column] = correct_column(fields[target_column])
                dest.write('%s\n' % ('\t'.join(fields)))

if __name__ == '__main__':
    correct_metadata(sys.stdin, sys.stdout)

