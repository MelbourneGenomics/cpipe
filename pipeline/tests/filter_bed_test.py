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

import unittest
import imp
import os
import random
import re
import sys
import io

sys.path.append('../scripts/')
import filter_bed

class FilterBedTest(unittest.TestCase):

    def test_negative(self):
        source = ['chr1\t100\t200\tA\t1\t0\n', 'chr1\t150\t100\tA\t2\t0\n', 'chr1\t200\t300\tA\t3\t10\n']
        target = io.StringIO()
        filter_bed.filter_bed(source, target)
        lines = target.getvalue().split('\n')
        assert lines[0] == 'chr1\t100\t200\tA\t1\t0'
        assert lines[1] == 'chr1\t200\t300\tA\t3\t10'
        assert len(lines) == 3 # one empty

    def test_exclude(self):
        source = ['chr1\t100\t200\tA\t1\t0\n', 'chr1\t150\t100\tB\t2\t0\n', 'chr1\t200\t300\tC\t3\t10\n']
        target = io.StringIO()
        exclude = ['A\n', 'X\n']
        filter_bed.filter_bed(source, target, exclude)
        lines = target.getvalue().split('\n')
        assert lines[0] == 'chr1\t200\t300\tC\t3\t10'
        assert len(lines) == 2 # one empty

    def test_include(self):
        source = ['chr1\t100\t200\tA\t1\t0\n', 'chr1\t150\t100\tB\t2\t0\n', 'chr1\t200\t300\tC\t3\t10\n']
        target = io.StringIO()
        include = ['C\n', 'X\n']
        filter_bed.filter_bed(source, target, None, include)
        lines = target.getvalue().split('\n')
        assert lines[0] == 'chr1\t200\t300\tC\t3\t10'
        assert len(lines) == 2 # one empty

    def test_no_genes(self):
        source = ['chr1\t100\t200\n', 'chr1\t200\t300\n']
        target = io.StringIO()
        filter_bed.filter_bed(source, target, None, None)
        lines = target.getvalue().split('\n')
        assert lines[0] == 'chr1\t100\t200'
        assert len(lines) == 3 # one empty
