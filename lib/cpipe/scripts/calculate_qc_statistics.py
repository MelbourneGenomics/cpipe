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
#
# Purpose:
#   Generate fragment statistics
# Usage:
#   samtools view file.bam | python calculate_qc_statistics.py > stats.out
###########################################################################
"""

import datetime
import sys

BASE_QUALITY_THRESHOLD = 30


def write_log(log, msg):
    """
        write a date stamped message to log
    """
    now = datetime.datetime.now().strftime('%y%m%d-%H%M%S')
    if log is not None:
        log.write('%s: %s\n' % (now, msg))


def calculate_statistics(sam, log):
    """
        read lines from sam and calculate mean and sd
    """
    write_log(log, 'started reading sam...')
    fragment_count = 0
    fragment_mean = 0.0
    fragment_std = 0.0
    read_count = 0
    read_mean = 0.0
    read_std = 0.0
    base_count = 0
    pass_count = 0
    for idx, line in enumerate(sam):
        fields = line.strip('\n').split('\t')
        flag = int(fields[1])
        # read
        incoming = float(len(fields[9]))
        read_count += 1
        if read_count == 1:
            read_mean = incoming
            read_std = 0.0
        else:
            old_mean = read_mean  # old m
            read_mean = old_mean + (incoming - old_mean) / read_count
            read_std += (incoming - old_mean) * (incoming - read_mean)

        # fragment
        if flag & 0x01 > 0 and flag & 0x02 > 0 and flag & 0x04 == 0:  # multiple segments, properly aligned, segment mapped
            incoming = float(fields[8])
            if incoming > 0:  # only +ve distances
                fragment_count += 1
                # calculate running sd calculation
                if fragment_count == 1:
                    fragment_mean = incoming
                    fragment_std = 0.0
                else:
                    old_mean = fragment_mean  # old mean
                    fragment_mean = old_mean + (incoming - old_mean) / fragment_count
                    fragment_std += (incoming - old_mean) * (incoming - fragment_mean)

        # quality (% of bases above given quality)
        incoming = [ord(c) - 33 for c in fields[10]]
        base_count += len(incoming)
        pass_count += sum([1 if c >= BASE_QUALITY_THRESHOLD else 0 for c in incoming])

        if idx % 100000 == 0:
            if fragment_count > 1:
                write_log(log, 'processed {0} pairs out of {1} lines. {2} ({3})...'.format(fragment_count, idx,
                                                                                           fragment_mean,
                                                                                           fragment_std / (
                                                                                               fragment_count - 1)))
            else:
                write_log(log, 'processed {0} lines...'.format(idx))

    # done

    # reads
    if read_count > 1:
        result = {'read_count': read_count, 'read_mean': read_mean, 'read_sd': (read_std / (read_count - 1)) ** 0.5}
    else:
        write_log(log, 'ERROR: no reads found')
        result = {'read_count': read_count, 'read_mean': 0.0, 'read_sd': 0.0}

    # fragments
    if fragment_count > 1:
        result.update({
            'fragment_count': fragment_count, 'fragment_mean': fragment_mean,
            'fragment_sd': (fragment_std / (fragment_count - 1)) ** 0.5
        })
    else:
        write_log(log, 'ERROR: no positive mate distances found')
        result.update({'fragment_count': fragment_count, 'fragment_mean': 0.0, 'fragment_sd': 0.0})

    # base quality
    result['base_count'] = base_count
    result['base_pass'] = pass_count

    return result


def main(sam=sys.stdin, out=sys.stdout, log=sys.stderr):
    """
        calculate fragment stats and write in tab separated format
    """
    stats = calculate_statistics(sam, log)
    out.write('fragment_count\t{0}\n'.format(stats['fragment_count']))
    out.write('fragment_mean\t{0}\n'.format(stats['fragment_mean']))
    out.write('fragment_sd\t{0}\n'.format(stats['fragment_sd']))
    out.write('read_count\t{0}\n'.format(stats['read_count']))
    out.write('read_mean\t{0}\n'.format(stats['read_mean']))
    out.write('read_sd\t{0}\n'.format(stats['read_sd']))
    out.write('base_count\t{0}\n'.format(stats['base_count']))
    out.write('base_pass\t{0}\n'.format(stats['base_pass']))


if __name__ == '__main__':
    main()
