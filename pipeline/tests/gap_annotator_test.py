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
import gap_annotator

class GapAnnotatorTest(unittest.TestCase):

    def test_simple_no_transcript(self):
        cov = ['chr1\t100\t200\tA\t1\t0', 'chr1\t100\t200\tA\t2\t0', 'chr1\t100\t200\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', '703\tNM_033083\tchr3\t+\t15469063\t15484120\t15469286\t15480662\t6\t15469063,15471419,15473593,15475854,15477848,15480615,\t15469389,15471514,15473730,15476045,15478082,15484120,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,']
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 0, target, ds, log)
        lines = target.getvalue().split('\n')
        assert lines[1] == 'chr1,A,100,101,0,0,0.0,2,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A'
        assert len(lines) == 3
        assert lines[2] == ''

    def test_simple_no_overlap(self):
        cov = ['chr1\t100\t200\tA\t1\t0', 'chr1\t100\t200\tA\t2\t0', 'chr1\t100\t200\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', '703\tNM_033083\tchr1\t+\t15469063\t15484120\t15469286\t15480662\t6\t15469063,15471419,15473593,15475854,15477848,15480615,\t15469389,15471514,15473730,15476045,15478082,15484120,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,']
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 0, target, ds, log)
        lines = target.getvalue().split('\n')
        assert lines[1] == 'chr1,A,100,101,0,0,0.0,2,NM_033083,+,15469184,N/A,N/A,N/A,N/A,1,N/A'
        assert len(lines) == 3
        assert lines[2] == ''

    def test_simple_overlap(self):
        cov = ['chr1\t15471419\t15471614\tA\t1\t0', 'chr1\t15471419\t15471614\tA\t2\t0', 'chr1\t15471419\t15471614\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', '703\tNM_033083\tchr1\t+\t15469063\t15484120\t15469286\t15480662\t6\t15469063,15471419,15473593,15475854,15477848,15480615,\t15469389,15471514,15473730,15476045,15478082,15484120,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,']
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 0, target, ds, log)
        lines = target.getvalue().split('\n')
        # Chr,Gene,Start,End,Min Cov,Max Cov,Median Cov,Width,Tx Name,Strand,CDS Distance,CDS Overlap Start,CDS Overlap End,AA Overlap Start,AA Overlap End,Exon Number,Exon Rank
        assert lines[1] == 'chr1,A,15471419,15471420,0,0,0.0,2,NM_033083,+,0,1,2,1,1,2,2'
        assert len(lines) == 3
        assert lines[2] == ''
