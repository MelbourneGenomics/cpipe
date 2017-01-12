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

import filter_transcripts

class FilterTranscriptsTest(unittest.TestCase):

    def test_simple_xm_only(self):
        '''
            test that XM on its own is not filtered
        '''
        unfiltered = ['CHROM\tPOS\tREF\tALT\tFeature\n', 'chr1\t100\tA\tA\tXM_1\n']
        log = io.StringIO()
        target = io.StringIO()
        filter_transcripts.filter_tsv(unfiltered, target, log)
        lines = target.getvalue().split('\n')
        assert lines[1] == 'chr1\t100\tA\tA\tXM_1'
        assert len(lines) == 3
        assert lines[2] == ''

    def test_simple_nm_only(self):
        '''
            test that NM on its own is not filtered
        '''
        unfiltered = ['CHROM\tPOS\tREF\tALT\tFeature\n', 'chr1\t100\tA\tA\tNM_1\n']
        log = io.StringIO()
        target = io.StringIO()
        filter_transcripts.filter_tsv(unfiltered, target, log)
        lines = target.getvalue().split('\n')
        assert lines[1] == 'chr1\t100\tA\tA\tNM_1'
        assert len(lines) == 3
        assert lines[2] == ''

    def test_simple_both(self):
        '''
            test that XM is filtered if an NM is present
        '''
        unfiltered = ['CHROM\tPOS\tREF\tALT\tFeature\n', 'chr1\t100\tA\tA\tNM_1\n', 'chr1\t100\tA\tA\tXM_1\n']
        log = io.StringIO()
        target = io.StringIO()
        filter_transcripts.filter_tsv(unfiltered, target, log)
        lines = target.getvalue().split('\n')
        assert lines[1] == 'chr1\t100\tA\tA\tNM_1'
        assert len(lines) == 3
        assert lines[2] == ''
