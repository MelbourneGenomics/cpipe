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
import filter_tsv

class FilterTSVTest(unittest.TestCase):

    def test_af(self):
        src = ['AF\tQUAL\ts1.DP\ts1.AD\n', '0.2\t10\t10\t3,7\n', '0.1\t10\t10\t3,7\n']
        log = StringIO.StringIO()
        target = StringIO.StringIO()
        filter_tsv.filter_tsv(src, target, log, af_min=0.15, qual_min=5, dp_min=5, ad_min=2, reverse=False)
        lines = target.getvalue().split('\n')
        assert lines[1] == '0.2\t10\t10\t3,7'
        assert lines[2] == ''
        assert len(lines) == 3

    def test_dp(self):
        src = ['AF\tQUAL\ts1.DP\ts1.AD\n', '0.2\t10\t1\t3,7\n', '0.2\t10\t10\t3,7\n']
        log = StringIO.StringIO()
        target = StringIO.StringIO()
        filter_tsv.filter_tsv(src, target, log, af_min=0.15, qual_min=5, dp_min=5, ad_min=2, reverse=False)
        lines = target.getvalue().split('\n')
        assert lines[1] == '0.2\t10\t10\t3,7'
        assert lines[2] == ''
        assert len(lines) == 3

    def test_ad_trio(self):
        src = ['AF\tQUAL\ts1.DP\ts1.AD\ts2.DP\ts2.AD\n', '0.2\t10\t1\t3,7\t10\t2,8\n', '0.2\t10\t7\t3,1\t4\t2,8\n']
        log = StringIO.StringIO()
        target = StringIO.StringIO()
        filter_tsv.filter_tsv(src, target, log, af_min=0.15, qual_min=5, dp_min=5, ad_min=2, reverse=False)
        lines = target.getvalue().split('\n')
        assert lines[1] == '0.2\t10\t1\t3,7\t10\t2,8'
        assert lines[2] == ''
        assert len(lines) == 3

    def test_dp_trio(self):
        src = ['AF\tQUAL\ts1.DP\ts1.AD\ts2.DP\ts2.AD\n', '0.2\t10\t1\t3,7\t10\t2,8\n', '0.2\t10\t5\t3,7\t4\t2,8\n']
        log = StringIO.StringIO()
        target = StringIO.StringIO()
        filter_tsv.filter_tsv(src, target, log, af_min=0.15, qual_min=5, dp_min=5, ad_min=2, reverse=False)
        lines = target.getvalue().split('\n')
        assert lines[1] == '0.2\t10\t1\t3,7\t10\t2,8'
        assert lines[2] == ''
        assert len(lines) == 3

    def test_dp_trio_proband_ok(self):
        src = ['AF\tQUAL\ts1.DP\ts1.AD\ts2.DP\ts2.AD\n', '0.2\t10\t5\t3,7\t4\t2,8\n', '0.2\t10\t4\t3,7\t5\t2,8\n']
        log = StringIO.StringIO()
        target = StringIO.StringIO()
        filter_tsv.filter_tsv(src, target, log, af_min=0.15, qual_min=5, dp_min=5, ad_min=2, reverse=False, proband='s1')
        lines = target.getvalue().split('\n')
        assert lines[1] == '0.2\t10\t5\t3,7\t4\t2,8'
        assert lines[2] == ''
        assert len(lines) == 3

    def test_ad(self):
        src = ['AF\tQUAL\ts1.DP\ts1.AD\n', '0.2\t10\t1\t3,1\n', '0.2\t10\t10\t3,7\n']
        log = StringIO.StringIO()
        target = StringIO.StringIO()
        filter_tsv.filter_tsv(src, target, log, af_min=0.15, qual_min=5, dp_min=5, ad_min=2, reverse=False)
        lines = target.getvalue().split('\n')
        assert lines[1] == '0.2\t10\t10\t3,7'
        assert lines[2] == ''
        assert len(lines) == 3
