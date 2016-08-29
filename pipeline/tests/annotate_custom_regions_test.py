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
import annotate_custom_regions

class AnnotateCustomRegionsTest(unittest.TestCase):

    def test_simple(self):
        regions = ['chr1\t100\t200\tA', 'chr1\t150\t250\tB', 'chr1\t200\t300\tC']
        log = StringIO.StringIO()
        source = ['CHROM\tPOS\tX', 'chr1\t100\tZ', 'chr1\t200\tY', 'chr1\t500\tW']
        target = StringIO.StringIO()
        annotate_custom_regions.annotate(source=source, target=target, regions=annotate_custom_regions.build_regions(regions, log), log=log)
        lines = target.getvalue().split('\n')
        assert lines[0] == 'CHROM\tPOS\tX\tCPIPE_BED'
        assert lines[1] == 'chr1\t100\tZ\tA'
        assert lines[2] == 'chr1\t200\tY\tB;C'
        assert lines[3] == 'chr1\t500\tW\t'

