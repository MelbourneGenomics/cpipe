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
        assert lines[1] == 'chr1,A,100,102,0,0,0.0,0.0,2,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A'
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
        assert lines[1] == 'chr1,A,100,102,0,0,0.0,0.0,2,NM_033083,+,15469185,N/A,N/A,N/A,N/A,N/A,N/A,1,1'
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
        # Chr,Gene,Start,End,Min Cov,Max Cov,Median Cov,Mean Cov,Width,Tx Name,Strand,CDS Distance,CDS Overlap Start,CDS Overlap End,AA Overlap Start,AA Overlap End,Exon Number,Exon Rank
        assert lines[1] == 'chr1,A,15471419,15471421,0,0,0.0,0.0,2,NM_033083,+,0,1,2,104,105,1,1,2,2'
        assert len(lines) == 3
        assert lines[2] == ''

    def test_simple_mean(self):
        cov = ['chr1\t15471419\t15471614\tA\t1\t0', 'chr1\t15471419\t15471614\tA\t2\t0', 'chr1\t15471419\t15471614\tA\t3\t1', 'chr1\t15471419\t15471614\tA\t4\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', '703\tNM_033083\tchr1\t+\t15469063\t15484120\t15469286\t15480662\t6\t15469063,15471419,15473593,15475854,15477848,15480615,\t15469389,15471514,15473730,15476045,15478082,15484120,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,']
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 1, target, ds, log)
        lines = target.getvalue().split('\n')
        # Chr,Gene,Start,End,Min Cov,Max Cov,Median Cov,Mean Cov,Width,Tx Name,Strand,CDS Distance,CDS Overlap Start,CDS Overlap End,AA Overlap Start,AA Overlap End,Exon Number,Exon Rank
        assert lines[1] == 'chr1,A,15471419,15471422,0,1,0,0.3,3,NM_033083,+,0,1,3,104,106,1,1,2,2'
        assert len(lines) == 3
        assert lines[2] == ''

    def test_multiple_overlap(self):
        cov = ['chr1\t15471419\t15471614\tA\t1\t0', 'chr1\t15471419\t15471614\tA\t2\t0', 'chr1\t15471419\t15471614\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', \
            '703\tNM_033083\tchr1\t+\t15469063\t15484120\t15469286\t15480662\t6\t15469063,15471419,15473593,15475854,15477848,15480615,\t15469389,15471514,15473730,15476045,15478082,15484120,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,', \
            '703\tdummy\tchr1\t+\t15469063\t15484120\t15469286\t15480662\t2\t15469063,15471412,\t15469389,15471516,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,' ]
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 0, target, ds, log)
        lines = target.getvalue().split('\n')
        # Chr,Gene,Start,End,Min Cov,Max Cov,Median Cov,Mean Cov,Width,Tx Name,Strand,CDS Distance,CDS Overlap Start,CDS Overlap End,AA Overlap Start,AA Overlap End,Exon Number,Exon Rank
        # can arrive in either order
        assert 'chr1,A,15471419,15471421,0,0,0.0,0.0,2,NM_033083,+,0,1,2,104,105,1,1,2,2' in lines
        assert 'chr1,A,15471419,15471421,0,0,0.0,0.0,2,dummy,+,0,8,9,111,112,3,3,2,2' in lines
        assert len(lines) == 4
        assert lines[3] == ''

    def test_exclude(self):
        cov = ['chr1\t15471419\t15471614\tA\t1\t0', 'chr1\t15471419\t15471614\tA\t2\t0', 'chr1\t15471419\t15471614\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', \
            '703\tNM_033083\tchr1\t+\t15469063\t15484120\t15469286\t15480662\t6\t15469063,15471419,15473593,15475854,15477848,15480615,\t15469389,15471514,15473730,15476045,15478082,15484120,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,', \
            '703\tdummy\tchr1\t+\t15469063\t15484120\t15469286\t15480662\t2\t15469063,15471412,\t15469389,15471516,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,' ]
        exclude = ['NM_033083',]
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log, exclude)
        gap_annotator.find_gaps(cov, 0, 0, target, ds, log)
        lines = target.getvalue().split('\n')
        # Chr,Gene,Start,End,Min Cov,Max Cov,Median Cov,Mean Cov,Width,Tx Name,Strand,CDS Distance,CDS Overlap Start,CDS Overlap End,AA Overlap Start,AA Overlap End,Exon Number,Exon Rank
        assert lines[1] == 'chr1,A,15471419,15471421,0,0,0.0,0.0,2,dummy,+,0,8,9,111,112,3,3,2,2'
        assert len(lines) == 3
        assert lines[2] == ''

    def test_distance(self):
        cov = ['chr1\t100\t200\tA\t1\t0', 'chr1\t100\t200\tA\t2\t0', 'chr1\t100\t200\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', \
            '703\tNM_033083\tchr1\t+\t100\t1000\t105\t920\t2\t103,500,\t180,550,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,' ]
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 0, target, ds, log)
        lines = target.getvalue().split('\n')
        assert lines[1] == 'chr1,A,100,102,0,0,0.0,0.0,2,NM_033083,+,4,N/A,N/A,N/A,N/A,N/A,N/A,1,1'
 
    def test_neg_distance(self):
        cov = ['chr1\t925\t1200\tA\t1\t0', 'chr1\t925\t1200\tA\t2\t0', 'chr1\t925\t1200\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', \
            '703\tNM_033083\tchr1\t+\t100\t1000\t105\t920\t2\t103,500,\t180,930,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,' ]
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 0, target, ds, log)
        lines = target.getvalue().split('\n')
        assert lines[1] == 'chr1,A,925,927,0,0,0.0,0.0,2,NM_033083,+,-6,N/A,N/A,N/A,N/A,N/A,N/A,2,2'
 
    def test_neg_strand(self):
        cov = ['chr1\t925\t1200\tA\t1\t0', 'chr1\t925\t1200\tA\t2\t0', 'chr1\t925\t1200\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', \
            '703\tNM_033083\tchr1\t-\t100\t1000\t105\t920\t2\t103,500,\t180,930,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,' ]
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 0, target, ds, log)
        lines = target.getvalue().split('\n')
        assert lines[1] == 'chr1,A,925,927,0,0,0.0,0.0,2,NM_033083,-,6,N/A,N/A,N/A,N/A,N/A,N/A,2,1'

    def test_cds_overlap(self):
        cov = ['chr1\t120\t1200\tA\t1\t0', 'chr1\t120\t1200\tA\t2\t0', 'chr1\t120\t1200\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', \
            '703\tNM_033083\tchr1\t-\t100\t1000\t105\t920\t2\t103,500,\t180,930,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,' ]
        target = StringIO.StringIO()
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 1, target, ds, log)
        lines = target.getvalue().split('\n')
        # Chr,Gene,Start,End,Min Cov,Max Cov,Median Cov,Mean Cov,Width,Tx Name,Strand,CDS Distance,CDS Overlap Start,CDS Overlap End,CDS Segment Start, CDS Segment End, AA Overlap Start,AA Overlap End,Exon Number,Exon Rank
        assert lines[1] == 'chr1,A,120,122,0,0,0.0,0.0,2,NM_033083,-,0,59,60,479,480,20,20,1,2' #chr1,A,120,121,0,1,0,0.3,3,NM_033083,+,0,1,3,436,437,1,1,2,2'
        assert len(lines) == 3
        assert lines[2] == ''

    def find_value(self, lines, colname):
        result = []
        index = lines[0].split(',').index(colname)
        for line in lines[1:]:
            fields = line.split(',')
            if len(fields) < index:
                result.append('')
            else:
                result.append(line.split(',')[index])
        return result

    def test_bed_overlap(self):
        cov = ['chr1\t120\t1200\tA\t1\t0', 'chr1\t120\t1200\tA\t2\t0', 'chr1\t120\t1200\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', \
            '703\tNM_033083\tchr1\t-\t100\t1000\t105\t920\t2\t103,500,\t180,930,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,' ]
        bed = ['chr1\t110\t130\tbed1']
        target = StringIO.StringIO()
        beds = {}
        beds['b1'] = gap_annotator.init_bed(bed, 'b1', log)
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 1, target, ds, log, beds=beds)
        lines = target.getvalue().split('\n')
        assert self.find_value(lines, 'b1')[0] == 'bed1'
        
    def test_bed_multi_overlap(self):
        cov = ['chr1\t120\t1200\tA\t1\t0', 'chr1\t120\t1200\tA\t2\t0', 'chr1\t120\t1200\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', \
            '703\tNM_033083\tchr1\t-\t100\t1000\t105\t920\t2\t103,500,\t180,930,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,' ]
        b1 = ['chr1\t110\t130\tbed1', 'chr1\t115\t125\tbed2']
        b2 = ['chr1\t110\t130\tb3']
        target = StringIO.StringIO()
        beds = {}
        beds['b1'] = gap_annotator.init_bed(b1, 'b1', log)
        beds['b2'] = gap_annotator.init_bed(b2, 'b2', log)
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 1, target, ds, log, beds=beds)
        lines = target.getvalue().split('\n')
        value = self.find_value(lines, 'b1')[0]
        assert 'bed1' in value
        assert 'bed2' in value
        assert self.find_value(lines, 'b2')[0] == 'b3'

    def test_bed_no_overlap(self):
        cov = ['chr1\t120\t1200\tA\t1\t0', 'chr1\t120\t1200\tA\t2\t0', 'chr1\t120\t1200\tA\t3\t10']
        log = StringIO.StringIO()
        data = ['bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames', \
            '703\tNM_033083\tchr1\t-\t100\t1000\t105\t920\t2\t103,500,\t180,930,\t0\tEAF1\tcmpl\tcmpl\t0,1,0,2,1,1,' ]
        bed = ['chr1\t210\t230\tbed1']
        target = StringIO.StringIO()
        beds = {}
        beds['b1'] = gap_annotator.init_bed(bed, 'b1', log)
        ds = gap_annotator.init_db(data, log)
        gap_annotator.find_gaps(cov, 0, 1, target, ds, log, beds=beds)
        lines = target.getvalue().split('\n')
        assert self.find_value(lines, 'b1')[0] == ''
