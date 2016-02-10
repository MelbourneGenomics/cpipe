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
import calculate_fragment_statistics

class CalculateFragmentTest(unittest.TestCase):

    def test_calc(self):
        bam = ['0\t3\t2\t3\t4\t5\t6\t7\t100\t9\t10\n', '0\t3\t2\t3\t4\t5\t6\t7\t-100\t9\t10\n', '0\t3\t2\t3\t4\t5\t6\t7\t200\t9\t10\n', '0\t3\t2\t3\t4\t5\t6\t7\t-200\t9\t10\n']
        log = StringIO.StringIO()
        result = calculate_fragment_statistics.calculate_statistics(bam, log)
        assert result['count'] == 2
        assert result['mean'] == 150.0
        assert int(result['sd'] * 100) == 7071

    def test_out(self):
        bam = ['0\t3\t2\t3\t4\t5\t6\t7\t100\t9\t10\n', '0\t3\t2\t3\t4\t5\t6\t7\t200\t9\t10\n', '0\t3\t2\t3\t4\t5\t6\t7\t300\t9\t10\n', '0\t3\t2\t3\t4\t5\t6\t7\t350\t9\t10\n']
        log = StringIO.StringIO()
        out = StringIO.StringIO()
        calculate_fragment_statistics.main(bam, out, log)
        lines = out.getvalue().split('\n')
        assert lines[0] == 'count\t4'
        assert lines[1] == 'mean\t237.5'
        assert lines[2] == 'sd\t110.86778913'

    def test_filtered(self):
        # should only include lines 1 and 4
        bam = ['0\t3\t2\t3\t4\t5\t6\t7\t100\t9\t10\n', '0\t7\t2\t3\t4\t5\t6\t7\t500\t9\t10\n', '0\t2\t2\t3\t4\t5\t6\t7\t1000\t9\t10\n', '0\t3\t2\t3\t4\t5\t6\t7\t200\t9\t10\n']
        log = StringIO.StringIO()
        result = calculate_fragment_statistics.calculate_statistics(bam, log)
        assert result['count'] == 2
        assert result['mean'] == 150.0
        assert int(result['sd'] * 100) == 7071
