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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.    If not, see <http:#www.gnu.org/licenses/>.
#
###########################################################################
# This script takes a sample metadata file as input (stdin) and
# prints out the field values in a human readable form
###########################################################################
'''

import sys

IS_NUMERIC = set(['dna concentration', 'dna quantity', 'dna quality', 'mean coverage'])
IS_ENUMERATION = {'sex': set(['Male', 'Female', 'Unknown', 'other']), 'sample type': set(['Normal', 'Tumour']), 'consanguinity': set(['No', 'Yes', 'Suspected', 'Unknown']), 'ethnicity': set(['Unknown', 'European', 'African', 'Asian'])}
IS_DATE = set(['dna_date', 'capture_date', 'sequencing_date'])

def is_valid_numeric(field):
    '''
        is numeric field valid
    '''
    if len(field) > 0:
        try:
            float(field)
            return True
        except ValueError:
            return False
    return True

def is_valid_enumeration(field, allowed):
    '''
        is enumerated field valid
    '''
    return len(field) == 0 or field in allowed

def is_valid_date(field):
    '''
        is date field valid
    '''
    return len(field) == 0 or len(field) == 8 and field.isdigit()

def validate(sample_fh, out, err):
    '''
        validate incoming fh
    '''
    headers = sample_fh.readline().strip().split('\t')
    idx = -1
    warnings = []
    for idx, line in enumerate(sample_fh):
        out.write("===== Sample {0} =====\n".format(idx))
        fields = line.strip('\n').split('\t')
        for jdx, field in enumerate(fields):
            out.write("{0:>24}: {1}\n".format(headers[jdx], field))
            if field.startswith(' '):
                warnings.append('Sample {0} field {1} (column {2}) contains leading whitespace'.format(idx, headers[jdx], jdx))
            if field.endswith(' '):
                warnings.append('Sample {0} field {1} (column {2}) contains trailing whitespace'.format(idx, headers[jdx], jdx))
            if headers[jdx].lower() in IS_NUMERIC and not is_valid_numeric(field):
                warnings.append('Sample {0} field {1} (column {2}) cannot be "{3}": must be empty or a number'.format(idx, headers[jdx], jdx, field))
            if headers[jdx].lower() in IS_ENUMERATION and not is_valid_enumeration(field, IS_ENUMERATION[headers[jdx].lower()]):
                warnings.append('Sample {0} field {1} (column {2}) cannot be "{3}": must be empty or one of: {4}'.format(idx, headers[jdx], jdx, field, ', '.join(IS_ENUMERATION[headers[jdx].lower()])))
            if headers[jdx].lower() in IS_DATE and not is_valid_date(field):
                warnings.append('Sample {0} field {1} (column {2}) cannot be "{3}": must be empty or date (yyyymmdd)'.format(idx, headers[jdx], jdx, field))

    if idx == -1:
        err.write("ERROR: file only contains one line. Are you using Windows style line feeds?\n")
    for warning in warnings:
        err.write("WARNING: {0}\n".format(warning))
    if len(warnings) == 0:
        err.write("No warnings\n")

if __name__ == '__main__':
    validate(sys.stdin, sys.stdout, sys.stderr)
