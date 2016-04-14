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
#
# Purpose:
#   Generate PED files from sample metadata file
# Usage:
#   See main() for arguments
#
###########################################################################
'''

import datetime
import sys

def write_log(fh, msg):
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    fh.write('{0}: {1}\n'.format(now, msg))

def generate_ped(sample_id, pedigree_info, gender, log):
    '''
        generate a ped file from the provided info
    '''
    lines = []
    lines.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format('#FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'))
    
    # extract pedigree info
    if '=' not in pedigree_info:
        write_log(log, 'ERROR: sample {0} has invalid pedigree {1}. Format should be fid=mid,pid'.format(sample_id, pedigree_info))
        return None

    # add members
    family_id, detail = pedigree_info.split('=')
    members = detail.split(',')
    parent = {}
    for member in members:
        lines.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(family_id, member, 0, 0, gender[member], 1))
        parent[gender[member]] = member

    # verify two parents of opposite gender
    if len(parent) != 2:
        write_log(log, 'ERROR: sample {0} has invalid pedigree {1}. There should be 2 parents but there were {2}'.format(sample_id, pedigree_info, len(parent)))
        return None
    if 'M' not in parent.values() or 'P' not in parent.values():
        write_log(log, 'ERROR: sample {0} has invalid pedigree {1}. Parents should have different genders: {2}'.format(sample_id, pedigree_info, parent))
        return None

    # add self
    lines.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(family_id, sample_id, parent[1], parent[2], gender[sample_id], 2))
    
    return lines

def generate_peds(metadata, log_fh):
    '''
        returns a dictionary of the form sample -> ped
    '''
    result = {}
    header = None
    sample_id = None
    pedigree_file = None
    gender_pos = None
    gender_map = {}
    count = 0
    found = 0
    lines = []
    write_log(log_fh, 'parsing sample metadata...')

    # get all genders
    for line in metadata:
        fields = line.strip('\n').split('\t')
        if header is None:
            header = fields
            sample_id = fields.index('Sample_ID')
            pedigree_file = fields.index('Pedigree_File')
            gender_pos = fields.index('Sex')
            continue
        gender = fields[gender_pos].strip() # male, female
        sample = fields[sample_id].strip()
        if gender.lower() == 'female':
            gender_map[sample] = 2 # female
        else:
            gender_map[sample] = 1 # male
        lines.append(line)

    # now build the pedigrees
    for count, line in enumerate(lines):
        fields = line.strip('\n').split('\t')
        sample = fields[sample_id].strip()
        pedigree = fields[pedigree_file].strip()
        if len(sample) > 0 and len(pedigree) > 0:
            result[sample] = generate_ped(sample, pedigree, gender_map, log_fh)
            found += 1

    write_log(log_fh, 'parsing sample metadata: processed {0} lines, found {1} families'.format(count, found))
    return result

def write_peds(prefix, data):
    '''
         write each item to the filesystem
    '''
    for key in data:
        with open('{0}{1}.ped'.format(prefix, key), 'w') as fh_out:
            if data[key] is not None:
                for line in data[key]:
                    fh_out.write(line)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Compare BAMs')
    parser.add_argument('--prefix', default='family_', help='prefix each file')
    args = parser.parse_args()

    peds = generate_peds(sys.stdin, sys.stderr)
    write_peds(args.prefix, peds)

if __name__ == '__main__':
    main()
