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

sys.path.append('../scripts/')
import extract_gatk_table_params

class ExtractGATKParamsTest(unittest.TestCase):

    def test_get_info(self):
        vcf = ['##INFO=<ID=ANN,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT">']
        out = io.StringIO()
        params = extract_gatk_table_params.extract_parameters(vcf, out)
        assert out.getvalue() == '-F ANN'

    def test_get_format(self):
        vcf = ['##FORMAT=<ID=ANN,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT">']
        out = io.StringIO()
        params = extract_gatk_table_params.extract_parameters(vcf, out)
        assert out.getvalue() == '-GF ANN'

