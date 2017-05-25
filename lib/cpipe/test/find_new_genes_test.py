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
import io

from cpipe.scripts import find_new_genes

class FindNewGenesTest(unittest.TestCase):

    def test_find(self):
      reference_genes = io.StringIO( 'c1\t1\t2\tabc\nc1\t3\t4\tdef\nc1\t5\t6\tghi' )
      excluded_genes = io.StringIO( 'ghi' )
      sample_lines = ['Sample_ID\tCohort\tPrioritised_Genes', '123\tCS\t4:def,ghi,jkl']
      log = io.StringIO()
      result = find_new_genes.generate_new_genes( sample_lines, log, reference_genes, excluded_genes, 'ref', 'exc' )
      assert list( result['CS']['notfound'] ) == [ 'JKL' ]
      assert list( result['CS']['add'] ) == [ 'DEF' ]
      assert list( result['CS']['addonce.123'] ) == [ 'GHI' ]

