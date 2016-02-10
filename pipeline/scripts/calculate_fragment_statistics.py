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
#   Generate fragment statistics
# Usage:
#   samtools view file.bam | python calculate_fragment_statistics.py > stats.out
###########################################################################
'''

import datetime
import sys

def write_log(log, msg):
    '''
        write a date stamped message to log
    '''
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    if log is not None:
        log.write('%s: %s\n' % (now, msg))

def calculate_statistics(sam, log):
    '''
        read lines from sam and calculate mean and sd
    '''
    write_log(log, 'started reading sam...')
    count = 0
    mean = 0.0
    std = 0.0
    for idx, line in enumerate(sam):
        fields = line.strip('\n').split('\t')
        flag = int(fields[1])
        if flag & 0x01 > 0 and flag & 0x02 > 0 and flag & 0x04 == 0:
            incoming = float(fields[8])
            if incoming > 0:
                count += 1
                # calculate running sd calculation
                if count == 1:
                    mean = incoming
                    std = 0.0
                else:
                    old_mean = mean # old m
                    mean = old_mean + (incoming - old_mean) / count
                    std = std + (incoming - old_mean) * (incoming - mean)
    
        if idx % 100000 == 0:
            if count > 1:
                write_log(log, 'processed {0} out of {1} lines. {2} ({3})...'.format(count, idx, mean, std / (count - 1)))
            else:
                write_log(log, 'processed {0} lines...'.format(idx))

    # done
    if count > 1:
        return {'count': count, 'mean': mean, 'sd': (std / (count - 1)) ** 0.5}
    else:
        write_log(log, 'ERROR: no positive mate distances found')
        return {'count': count, 'mean': 0.0, 'sd': 0.0}

def main(sam, out, log):
    '''
        calculate fragment stats and write in tab separated format
    '''
    stats = calculate_statistics(sam, log)
    out.write('count\t{0}\n'.format(stats['count']))
    out.write('mean\t{0}\n'.format(stats['mean']))
    out.write('sd\t{0}\n'.format(stats['sd']))

if __name__ == '__main__':
    main(sys.stdin, sys.stdout, sys.stderr)

