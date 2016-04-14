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
import convert_to_lovd

class ConvertToLovdTest(unittest.TestCase):

    def test_get_ann(self):
        vcf = ['##INFO=<ID=ANN,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT">']
        log = StringIO.StringIO()
        ann = convert_to_lovd.get_ann_fields(vcf, log)
        assert ann == ['Allele', 'Consequence', 'IMPACT']

    def test_table(self):
        table = []
        table.append('CHROM\tPOS\tID\tREF\tALT\tQUAL\tANN\t00NA12877.AB')
        table.append('chrM\t4746\t.\tA\tG\t6482.14\tG|intergenic_variant|MODIFIER\tNA')
        ann = ['Allele', 'Consequence', 'IMPACT']
        log = StringIO.StringIO()
        out = StringIO.StringIO()
        convert_to_lovd.process_table(table, out, ann, log)
        result = out.getvalue().split('\n')[1].split('\t') 
        assert len(result) == 10
        assert result == ['chrM', '4746', '.', 'A', 'G', '6482.14', 'G', 'intergenic_variant', 'MODIFIER', 'unknown']

    def test_no_ann(self):
        table = []
        table.append('CHROM\tPOS\tID\tREF\tALT\tQUAL\tANN\t00NA12877.AB')
        table.append('chrM\t4746\t.\tA\tG\t6482.14\t\tNA')
        ann = ['Allele', 'Consequence', 'IMPACT']
        log = StringIO.StringIO()
        out = StringIO.StringIO()
        convert_to_lovd.process_table(table, out, ann, log)
        result = out.getvalue().split('\n')[1].split('\t') 
        assert len(result) == 10
        assert result == ['chrM', '4746', '.', 'A', 'G', '6482.14', 'unknown', 'unknown', 'unknown', 'unknown']
