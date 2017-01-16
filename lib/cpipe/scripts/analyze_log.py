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
# Synopsis:
#   find stage times from log file
# Usage:
#   python analyze_log.py < logfile
#
###########################################################################
'''


import collections
import datetime
import sys

def analyze_log(log, result):
    '''
        looks for profile markers and generates a summary of time taken
    '''
    in_progress = collections.defaultdict(list)
    final = collections.defaultdict(list)
    current = 0
    max_in_progress = 0
    for line in log:
       # parse out enter and exit lines example:
       # 160505-033602: merge_bams: exit (020444901)
       if ": enter " in line:
           when, command, rest = line.split(': ')
           branch = rest.split('(')[-1].strip('()')
           identifier = '{0}/{1}'.format(command, branch)
           when = datetime.datetime(int(when[:2]), int(when[2:4]), int(when[4:6]), int(when[7:9]), int(when[9:11]), int(when[11:13]))
           in_progress[identifier].append(when)
           current += 1
           max_in_progress = max(max_in_progress, current)
       elif ": exit " in line:
           # find existing record
           when, command, rest = line.split(': ')
           branch = rest.split('(')[-1].strip('()')
           identifier = '{0}/{1}'.format(command, branch)
           if len(in_progress[identifier]) > 0:
               # determine time taken
               old_time = in_progress[identifier].pop(0)
               new_time = datetime.datetime(int(when[:2]), int(when[2:4]), int(when[4:6]), int(when[7:9]), int(when[9:11]), int(when[11:13]))
               final[command].append((new_time - old_time).total_seconds())
           else:
               sys.stderr.write('WARN: not found: ' + line)
           current -= 1
       else:
           pass
    # finished
    result.write('Stage\tMin\tMax\tAverage\tTotal\t%\n')
    total = sum([sum(final[key]) for key in list(final.keys())])
    for key in sorted(final.keys()):
        result.write('{0}\t{1}\t{2}\t{3:.0f}\t{4}\t{5:.0f}\n'.format(key, min(final[key]), max(final[key]), sum(final[key]) / len(final[key]), sum(final[key]), sum(final[key]) * 100 / total))
    result.write('Max concurrency: {0}\n'.format(max_in_progress))

if __name__ == '__main__':
    analyze_log(sys.stdin, sys.stdout)
