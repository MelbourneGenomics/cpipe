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
import genelist_to_bed

class GeneListToBedTest(unittest.TestCase):

    def test_exclude(self):
      bed_in = StringIO.StringIO( 'c1\t1\t2\tabc\nc1\t3\t4\tdef\nc1\t5\t6\tghi' )
      bed_out = StringIO.StringIO()
      log = StringIO.StringIO()
      genelist_to_bed.filter_bed( ['sample_genelist.txt'], bed_in, bed_out, log, exclude=['ghi'] )
      assert bed_out.getvalue() == 'c1\t1\t2\tabc\nc1\t3\t4\tdef\n'
      
    def test_notfound(self):
      bed_in = StringIO.StringIO( 'c1\t1\t2\tabc\nc1\t5\t6\tghi' )
      bed_out = StringIO.StringIO()
      log = StringIO.StringIO()
      genelist_to_bed.filter_bed( ['sample_genelist.txt'], bed_in, bed_out, log, exclude=['ghi'] )
      assert bed_out.getvalue() == 'c1\t1\t2\tabc\n'
      
