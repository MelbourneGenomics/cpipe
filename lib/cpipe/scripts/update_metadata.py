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
# Update fields in the metadata file
###########################################################################
"""

import argparse
import sys

def update_metadata(sample_in, sample_out, log, sample, name, value):
    """
        update sample with new value
    """
    if len(sample_in) == 0:
        log.write('ERROR: file is empty\n')
        return 1
    headers = sample_in[0].strip().split('\t')
    headers_map = {header.lower(): idx for idx, header in enumerate(headers)}
    sample_out.write(sample_in[0])

    if 'sample_id' not in headers_map:
        log.write('ERROR: source does not appear to be a valid metadata file\n')
        return 1

    if name.lower() not in headers_map:
        log.write('ERROR: {0} is not a valid field name. Must be one of: {1}\n'.format(name, ', '.join(list(headers_map.keys()))))
        return 1

    sample_found = False
    samples = []
    for line in sample_in[1:]:
        fields = line.rstrip('\n').split('\t')
        samples.append(fields[headers_map['sample_id']].strip())
        if samples[-1] == sample:
            sample_found = True
            fields[headers_map[name.lower()]] = value
            sample_found = True
            sample_out.write('%s\n' % '\t'.join(fields))
        else:
            sample_out.write(line)

    if not sample_found:
        log.write('ERROR: sample "{0}" not found in: {1}\n'.format(sample, ', '.join(samples)))

def main():
    """
        update metadata file from command line
    """
    parser = argparse.ArgumentParser(description='Update metadata')
    parser.add_argument('--sample_id', required=True, help='sample ID to update')
    parser.add_argument('--name', required=True, help='name of field to update')
    parser.add_argument('--value', required=True, help='new value for field')
    parser.add_argument('--target', required=True, help='filename')
    args = parser.parse_args()

    with open(args.target, 'r') as sample_in:
        lines = sample_in.readlines()
    with open(args.target, 'w') as sample_out:
        update_metadata(lines, sample_out, sys.stderr, args.sample_id, args.name, args.value)

if __name__ == '__main__':
    main()
