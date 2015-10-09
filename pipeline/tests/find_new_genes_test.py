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
import find_new_genes

class FindNewGenesTest(unittest.TestCase):

    def test_find(self):
      reference_genes = StringIO.StringIO( 'c1\t1\t2\tabc\nc1\t3\t4\tdef\nc1\t5\t6\tghi' )
      excluded_genes = StringIO.StringIO( 'ghi' )
      sample_lines = ['Cohort\tPrioritised_Genes', 'CS\t4:def,ghi,jkl']
      log = StringIO.StringIO()
      result = find_new_genes.generate_new_genes( sample_lines, log, reference_genes, excluded_genes )
      assert list( result['CS']['notfound'] ) == [ 'JKL' ]
      assert list( result['CS']['add'] ) == [ 'DEF' ]
      assert list( result['CS']['addonce'] ) == [ 'GHI' ]
