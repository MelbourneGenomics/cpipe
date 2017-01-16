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
'''

import argparse
import datetime
import os
import os.path
import sys

def run(cmd, log, dry=False):
    '''
        execute a command and log
    '''
    if dry:
        log.write('would execute: {0}\n'.format(cmd))
    else:
        log.write('executing: {0}\n'.format(cmd))
        os.system(cmd)

def mark_batch_finished(directory, log, read_only, move, dry=False):
    '''
        mark a batch read-only and move to completion
    '''
    # first mark as read only
    batch_dir = os.path.split(directory)[0] # e.g. ./batches/123
    if read_only:
        # we can't simply mark the directory read only because bpipe still needs to write things
        # e.g. find /hsm/VR0320/shared/pgeorgeson/cpipe-dev-deploy/cpipe-mgha/batches/trio_test -path "*/.bpipe" -prune -o -name commandlog.txt -prune -o -exec chmod -w "{}" \;
        #run('find "{0}" chmod -R -w "{0}"'.format(batch_dir), log, dry)
        run('find "{0}" -path "*/.bpipe" -prune -o -name commandlog.txt -prune -o -exec chmod -w "{{}}" \\;'.format(batch_dir), log, dry)
        

    batch_collection, batch_name = os.path.split(batch_dir) # ./batches, 123

    # now move the batch to a "completed" state
    if move:
        new_location = '{0}/complete.{1}.{2}'.format(batch_collection, batch_name, datetime.datetime.now().strftime('%y%m%d-%H%M%S'))
        run('mv {0} {1}'.format(batch_dir, new_location), log, dry)
        run('ln -s {1} {0}'.format(batch_dir, new_location), log, dry) # bpipe still needs to know about the old directory
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mark Cpipe batch as complete')
    parser.add_argument('--dry', required=False, default=False, action='store_true', help='dry run')
    parser.add_argument('--read_only', required=False, default=False, action='store_true', help='mark files read only')
    parser.add_argument('--move', required=False, default=False, action='store_true', help='move directory to completion')
    args = parser.parse_args()

    mark_batch_finished(os.getcwd(), sys.stderr, read_only=args.read_only, move=args.move, dry=args.dry)
