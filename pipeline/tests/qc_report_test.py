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
import qc_report

class QCReportTest(unittest.TestCase):

    def test_male(self):
        exome_cov = ['chr1\t100\t200\t1\t20', 'chrX\t100\t200\t1\t10', 'chrY\t100\t200\t1\t10']
        log = StringIO.StringIO()
        result = qc_report.calculate_karyotype(exome_cov, log)
        assert result['sex'] == 'MALE'
        assert result['x_mean_coverage'] == 10.
        assert result['y_mean_coverage'] == 10.
        assert result['autosome_mean_coverage'] == 20.

    def test_female(self):
        exome_cov = ['chr1\t100\t200\t1\t5', 'chrX\t100\t200\t1\t40', 'chrX\t400\t500\t1\t60']
        log = StringIO.StringIO()
        result = qc_report.calculate_karyotype(exome_cov, log)
        assert result['sex'] == 'FEMALE'
        assert result['x_mean_coverage'] == 50.
        assert result['y_mean_coverage'] == 0.
        assert result['autosome_mean_coverage'] == 5.

    def test_median(self):
        l = [ 5, 8, 2, 9, -3 ]
        assert qc_report.median(l) == 5
        l = [ -3, 8, 5, 8, 2, 9 ]
        assert qc_report.median(l) == 6.5

    def test_parse_metadata(self):
        m = [ 'Batch\tSample_ID\tDNA_Tube_ID\n', '026\t12345\t\n', '026\t11111\t22222\n' ]
        p = qc_report.parse_metadata(m, '12345')
        assert p['dna_tube_id'] == ''
        p = qc_report.parse_metadata(m, '11111')
        assert p['dna_tube_id'] == '22222'
        p = qc_report.parse_metadata(m, 'noway')
        assert p is None

    def test_calculate_mean_stats(self):
        log = StringIO.StringIO()
        mean = 20
        # expect first to include those in the range 4 -> 36 (80% either side)
        stats = [0, 10, 20, 100]
        result = qc_report.calculate_mean_stats(stats, mean, log)
        assert result == [50.0, 75.0, 75.0, 50.0, 25.0]

    def test_calculate_summary(self):
        cov = ['chr1\t100\t200\tA\t1\t5', 'chr1\t100\t200\tA\t2\t5', 'chr1\t100\t200\tA\t3\t15', 'chr1\t100\t200\tA\t4\t40']
        log = StringIO.StringIO()
        s = qc_report.calculate_summary(cov, 20, log)
        assert s['mean'] == 16.25
        assert s['median'] == 10
        assert s['genes']['A']['median'] == 10
        assert s['genes']['A']['ok'] == 25.0

#    def test_report(self):
#        cov = ['chr1\t100\t200\tA\t1\t5', 'chr1\t100\t200\tA\t2\t5', 'chr1\t100\t200\tA\t3\t15', 'chr1\t100\t200\tA\t4\t40']
#        log = StringIO.StringIO()
#        s = qc_report.calculate_summary(cov, 20, log)
#
#        exome_cov = ['chr1\t100\t200\t1\t5', 'chrX\t100\t200\t1\t40', 'chrX\t400\t500\t1\t60']
#        k = qc_report.calculate_karyotype(exome_cov, log)
#
#        m = [ 'Batch\tSample_ID\tDNA_Tube_ID\tSex\tMean_Coverage\n', '026\t12345\t\tFemale\t100.2\n', '026\t11111\t22222\t\n' ]
#        p = qc_report.parse_metadata(m, '12345')
#
#        threshold = 20
#
#        out = StringIO.StringIO()
#        qc_report.generate_report(s, k, p, threshold, out, log)
#        print out.getvalue()

    def test_parse_tsv(self):
        frags = ['mean\t123.4\n', 'sd\t567.9\n']
        res = qc_report.parse_tsv(frags)
        assert res['mean'] == '123.4'
        assert res['sd'] == '567.9'

    def test_parse_date(self):
        d = ''
        assert qc_report.parse_date(d) == 'N/A'
        d = '"20151230"'
        assert qc_report.parse_date(d) == '2015-12-30'
        d = '"20151230 20160112"'
        assert qc_report.parse_date(d) == '2015-12-30 2016-01-12'

    def test_build_metrics(self):
        p = ['## net.sf.picard.metrics.StringHeader', '## METRICS CLASS\tnet.sf.picard.sam.DuplicationMetrics', 'LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE', 'null\t79055\t45114896\t294963\t57840\t9210954\t2246706\t0.204628\t117222885', '## HISTOGRAM\tjava.lang.Double' ]
        o = ['12345']
        log = StringIO.StringIO()
        m = qc_report.build_metrics(p, o, log)
        assert m['read_pairs_examined'] == '45114896'

    def test_build_empty_categories(self):
        categories = ''
        prioritized = ''
        log = StringIO.StringIO()
        qc_report.build_categories(categories, prioritized, log)

    def test_grouping(self):
        assert qc_report.group_number(123) == '123'
        assert qc_report.group_number(1234) == '1,234'
        assert qc_report.group_number(1234567) == '1,234,567'
