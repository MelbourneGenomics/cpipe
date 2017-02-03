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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Cpipe.    If not, see <http:#www.gnu.org/licenses/>.
#
####################################################################################
#
# Purpose:
#     Look for excluded genes on a series of gene lists
# Usage:
#     validate_genelists.py --exclude exclusionfile genefile(s)...
####################################################################################
"""

import argparse
import sys

def find_excluded(exclude_fh, files, out):
    """
        given a file containing genes to exclude, read files and find genes matching exclusion list
    """
    excluded = set()
    for line in exclude_fh:
        if line.startswith('#'):
            continue
        excluded.add(line.split()[0].strip().upper())
    out.write('Testing {0} excluded genes\n\n'.format(len(excluded)))
    out.write('Cohort | Count | Genes \n')
    out.write('-------|-------|--------------------\n')

    total = set()
    for filename in files:
        included = set()
        with open(filename, 'r') as filehandle:
            for line in filehandle:
                if line.startswith('#'):
                    continue
                included.add(line.split()[0].strip().upper())
        result = included.intersection(excluded)
        out.write('{0} | {1} out of {2} | {3}\n'.format(filename, len(result), len(included), ' '.join(sorted(list(result)))))
        total = total.union(result)
    out.write('TOTAL | {0} excluded gene(s) found | {1}\n'.format(len(total), ' '.join(sorted(list(total)))))

def main():
    """
        execute from command line
    """
    parser = argparse.ArgumentParser(description='Validate gene list')
    parser.add_argument('--exclude', help='file containing genes to exclude')
    parser.add_argument('list', nargs='+', help='list of files to test')
    args = parser.parse_args()
    find_excluded(open(args.exclude, 'r'), args.list, sys.stdout)

if __name__ == '__main__':
    main()
