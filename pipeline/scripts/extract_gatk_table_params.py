#!/usr/bin/env python
'''
###########################################################################

 This file is part of Cpipe.

 Cpipe is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, under version 3 of the License, subject
 to additional terms compatible with the GNU General Public License version 3,
 specified in the LICENSE file that is part of the Cpipe distribution.

 Cpipe is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Cpipe.  If not, see <http:#www.gnu.org/licenses/>.

###########################################################################
Description:
    extract field names from a vcf for use in the GATK command VariantsToTable
Usage:
    python extract_gatk_table_params < vcf_file 
###########################################################################
'''

import re
import sys

def extract_parameters(vcf, target):
    '''
        extract field names from a vcf for use in the GATK command VariantsToTable
    '''
    result = []
    for line in vcf:
        if not line.startswith('#'):
            break # done
        # look for INFO field
        m = re.match(r"^##INFO=<ID=([^,]*)", line)
        if m is not None:
            result.append('-F {0}'.format(m.group(1)))
        m = re.match(r"^##FORMAT=<ID=([^,]*)", line)
        if m is not None:
            result.append('-GF {0}'.format(m.group(1)))
    target.write(' '.join(result))

if __name__ == '__main__':
    extract_parameters(sys.stdin, sys.stdout)
