#!/usr/bin/env python
"""
# This script has been deprecated and may be removed.
# Refer to combine_target_regions.py to generate bedfiles from gene lists
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
# Purpose:
#     Given some gene lists and reference bed file, generate a bed file with just those genes
# Usage:
#     genelist_to_bed genelist... < ref.bed > filtered.bed
#     optional arguments:
#     --exclude a file containing genes to exclude
####################################################################################
"""

import sys

def filter_bed(genelists, bed_in, bed_out, log, exclude=None):
    """ get the list of proposed genes """
    genes = set()
    for arg in genelists:
        for line in open(arg, 'r'):
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            genes.add(fields[0].upper())
    log.write('%i candidate genes added\n' % len(genes))

    # get the list of exclusions
    disallowed = set()
    if exclude is not None:
        for line in exclude:
            disallowed.add(line.strip().upper())
        log.write('%i excluded genes added\n' % len(disallowed))

    # filter the reference
    filtered = 0
    allowed = 0
    found = set()
    blocked = set()
    for line in bed_in:
        if line.startswith('#'):
            bed_out.write(line)
            allowed += 1
        else:
            fields = line.strip().split('\t')
            if len(fields) > 3:
                candidate = fields[3].upper()
                if candidate in genes:
                    if candidate in disallowed:
                        blocked.add(candidate)
                    else:
                        bed_out.write(line)
                        allowed += 1
                        found.add(candidate)
                else:
                    filtered += 1
            else:
                bed_out.write(line)
                allowed += 1

    log.write('%i lines written, %i lines filtered, %i out of %i candidate genes found\n' % (allowed, filtered, len(found), len(genes)))
    if len(found) != len(genes):
        log.write('Not found: %s\n' % ' '.join(sorted(list(genes.difference(found.union(blocked))))))
    if len(blocked) > 0:
        log.write('Excluded: %s\n' % ' '.join(sorted(list(blocked))))

def main():
    """
        run from command line
    """
    import argparse
    parser = argparse.ArgumentParser(description='Convert genelists to bed file')
    parser.add_argument('--exclude', dest='exclude', required=False, help='file containing genes to exclude')
    parser.add_argument('genelists', nargs='*', help='gene list files to include')
    args = parser.parse_args()
    if args.exclude:
        filter_bed(args.genelists, sys.stdin, sys.stdout, sys.stderr, exclude=open(args.exclude, 'r'))
    else:
        filter_bed(args.genelists, sys.stdin, sys.stdout, sys.stderr)

if __name__ == '__main__':
    main()
