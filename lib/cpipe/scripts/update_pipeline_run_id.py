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
'''

import argparse
import os.path
import random
import sys

def generate_new_id(filename):
    '''
        given a file, reads the current ID, appends to it, and writes it back to the same file.
        if the file doesn't exist, a random ID is generated.
        format of the ID is site_000000000
    '''
    current_id = get_current_id(filename)
    site, run = current_id.rsplit("_", 1)
    # increment run ID
    run = int(run) + 1

    new_id = '%s_%09i' % (site, run)

    handle = open(filename, 'w')
    handle.write(new_id)
    return new_id

def get_current_id(filename):
    '''
        given a file, reads the current ID. if the file doesn't exist, a random ID is generated.
        format of the ID is site_000000000
        @param: filename containing pipeline ID
        @returns: current ID
    '''
    # NOTE! we don't do any file locking.
    # parallel pipelines could potentially attempt to update the ID simulatenously, resulting in a non-unique ID
    if os.path.isfile(filename):
        handle = open(filename, 'r')
        current = handle.readline().strip()
        handle.close()
        if "_" in current:
            site, run = current.rsplit("_", 1)
            try:
                run = int(run)
            except ValueError: # rightmost section isn't an id after all
                site = current
                run = 1
            current_id = '%s_%09i' % (site, run)
        elif len(current) == 0: # empty file
            run = 0
            current_id = 'site%i_%09i' % (random.randint(0, 1e6), run)
        else: # no _
            current_id = '%s_%09i' % (current, 0)
    else: # no file
        run = 0
        current_id = 'site%i_%09i' % (random.randint(0, 1e6), run)

    return current_id

def write(src, target, new_id):
    '''
        reads lines from src and writes to target, appending the pipeline ID at the start as a new column of a tab separated file
    '''
    first = True
    for line in src:
        if first:
            target.write('Pipeline_Run_ID\t%s' % line)
            first = False
        else:
            target.write('%s\t%s' % (new_id, line))

def main():
    '''
        command line implementation
    '''
    parser = argparse.ArgumentParser(description='Generate sample metadata file with pipeline ID')
    parser.add_argument('--id', required=True, help='ID file to read/write')
    parser.add_argument('--increment', type=bool, required=False, default=False, help='Increment the pipeline ID')
    parser.add_argument('--parse', type=bool, required=False, default=False, help='Parse metadata file')
    args = parser.parse_args()
    if args.increment:
        new_pipeline_id = generate_new_id(args.id)
    else:
        new_pipeline_id = get_current_id(args.id)
    if args.parse:
        write(sys.stdin, sys.stdout, new_pipeline_id)
    else:
        print(new_pipeline_id)

if __name__ == "__main__":
    main()
