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

import generate_peds

class GeneratePedsTest(unittest.TestCase):
    def test_simple(self):
        metadata = ['Sample_ID\tPedigree_File\tSex\n', 'P\t\tMale\n', 'M\t\tFemale\n', 'C\tfid=P,M\tM\n']
        log = io.StringIO()
        result = generate_peds.generate_peds(metadata, log)
        #print 'log', log.getvalue()
        #print 'result', result
        assert result['C'] == ['#FID\tIID\tPID\tMID\tSex\tPhenotype\n', 'fid\tP\t0\t0\t1\t1\n', 'fid\tM\t0\t0\t2\t1\n', 'fid\tC\tP\tM\t1\t2\n']

    def test_unknown(self):
        metadata = ['Sample_ID\tPedigree_File\tSex\n', 'P\t\tUnknown\n', 'M\t\tUnknown\n', 'C\tfid=P,M\tM\n']
        log = io.StringIO()
        result = generate_peds.generate_peds(metadata, log)
        #print 'log', log.getvalue()
        #print 'result', result
        assert result['C'] == None
