#!/usr/bin/env python
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

import unittest
import imp
import os
import random
import re
import sys
import StringIO

sys.path.append('../scripts/')
import calculate_qc_statistics

class CalculateQCTest(unittest.TestCase):

    def test_calc(self):
        bam = ['0\t3\t2\t3\t4\t5\t6\t7\t100\tAAAA\tJJJJ\n', '0\t3\t2\t3\t4\t5\t6\t7\t-100\tAAAA\tJJJJ\n', '0\t3\t2\t3\t4\t5\t6\t7\t200\tAAAA\tJJJJ\n', '0\t3\t2\t3\t4\t5\t6\t7\t-200\tAAAA\tJJJJ\n']
        log = StringIO.StringIO()
        result = calculate_qc_statistics.calculate_statistics(bam, log)
        assert result['fragment_count'] == 2
        assert result['fragment_mean'] == 150.0
        assert int(result['fragment_sd'] * 100) == 7071
        assert result['read_count'] == 4
        assert result['read_mean'] == 4.0
        assert int(result['read_sd'] * 100) == 0
        assert result['base_count'] == 16
        assert result['base_pass'] == 16

    def test_out(self):
        bam = ['0\t3\t2\t3\t4\t5\t6\t7\t100\tA\t5\n', '0\t3\t2\t3\t4\t5\t6\t7\t200\tAB\t5I\n', '0\t3\t2\t3\t4\t5\t6\t7\t300\tABC\t555\n', '0\t3\t2\t3\t4\t5\t6\t7\t350\tABCD\t5II5\n']
        log = StringIO.StringIO()
        out = StringIO.StringIO()
        calculate_qc_statistics.main(bam, out, log)
        lines = out.getvalue().split('\n')
        assert lines[0] == 'fragment_count\t4'
        assert lines[1] == 'fragment_mean\t237.5'
        assert lines[2] == 'fragment_sd\t110.86778913'
        assert lines[3] == 'read_count\t4'
        assert lines[4] == 'read_mean\t2.5'
        assert lines[5] == 'read_sd\t1.29099444874'
        assert lines[6] == 'base_count\t10'
        assert lines[7] == 'base_pass\t3'

    def test_filtered(self):
        # should only include lines 1 and 4 in fragment count
        bam = ['0\t3\t2\t3\t4\t5\t6\t7\t100\t9\t10\n', '0\t7\t2\t3\t4\t5\t6\t7\t500\t9\t10\n', '0\t2\t2\t3\t4\t5\t6\t7\t1000\t9\t10\n', '0\t3\t2\t3\t4\t5\t6\t7\t200\t9\t10\n']
        log = StringIO.StringIO()
        result = calculate_qc_statistics.calculate_statistics(bam, log)
        assert result['fragment_count'] == 2
        assert result['fragment_mean'] == 150.0
        assert int(result['fragment_sd'] * 100) == 7071
