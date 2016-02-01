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
import StringIO

sys.path.append('../scripts/')
import gap_annotator

class GapAnnotatorTest(unittest.TestCase):

    def test_simple(self):
        cov = ['chr1\t100\t200\tA\t1\t0', 'chr1\t100\t200\tA\t2\t0', 'chr1\t100\t200\tA\t3\t10']
        log = StringIO.StringIO()
        target = StringIO.StringIO()
        gap_annotator.find_gaps(cov, 0, 0, target, {}, log)
        lines = target.getvalue().split('\n')
        assert lines[1] == 'chr1,A,100,101,0,0,0.0,2'
        assert len(lines) == 3
        assert lines[2] == ''
