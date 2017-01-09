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

import update_gene_lists

class UpdateGenesListsTest(unittest.TestCase):

    def test_find(self):
      if not os.path.exists('CS'):
        os.makedirs('CS')
      with open( 'CS/CS.genes.txt', 'w' ) as current:
        current.write( '#version\nabc\t1\ndef\t2\n' )
      with open( 'CS.add.genes.txt', 'w' ) as current:
        current.write( 'def\nghi\n' )
      log = io.StringIO()
      update_gene_lists.update_gene_lists( '.', '.', log )
      lines = open( 'CS/CS.genes.txt', 'r' ).readlines()
      assert lines[1:] == [ '#notes 1 gene(s) added: GHI\n', 'ABC\t1\n', 'DEF\t2\n', 'GHI\t1\n' ]
