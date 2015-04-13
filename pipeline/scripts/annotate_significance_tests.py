# vim: ts=4:expandtab:sw=4:cindent
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
import csv

#
# NOTE: to run these tests use  :
#
#    python annotate_significance_tests.py
#

from annotate_significance import Annovar 

class AnnotateSignificanceTest(unittest.TestCase):
    
    
    header = "Func,Gene,ExonicFunc,AAChange,Conserved,SegDup,ESP5400_ALL,1000g2010nov_ALL,dbSNP138,AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,Chr,Start,End,Ref,Obs,Otherinfo,Qual,Depth,Condel".split(",")
    
    # A prototype line that we use to create test data
    line = '"exonic","SCN5A","nonsynonymous SNV","NM_000335:c.G1339T:p.A447S","437;Name=lod=80",,,,,0.42,,,,,,,,,,,,chr3,38646399,38646399,C,A,"het","14.91","19","0.3"\n'

    def __init__(self,n):
        unittest.TestCase.__init__(self,n)
        Annovar.init_columns(self.header)
        variant = csv.reader([self.line], delimiter=",", quotechar='"').next()
        self.a = Annovar(variant)
        print "init ..."
     
    def testNovel(self):
        
        a = self.a
       
        # It has no MAF and no dbSNP ID => is novel
        assert a.is_novel()
        
        # But if we give it a dbSNP ID, then it is not novel
        a.set_value('dbSNP138','rs12345')
        assert not a.is_novel()
        a.set_value('dbSNP138','')
        
        assert a.is_novel()
        
        # Give it  1000 Genomes MAF and it should be not novel
        a.set_value('1000g2010nov_ALL','0.00001')
        assert not a.is_novel()

        
    def testPriorityIndex(self):
        
        a = self.a
        
        # Novel, missense, no Condel score but conserved by conserved region annotation: category 3
        a.set_value('Condel','')
        print a.priority()
        assert a.priority() == 3
        
        # Adding high Condel score makes it remain priority 3
        a.set_value('Condel','0.9')
        assert a.priority() == 3
        
        # Removing conserved annotation should not make a difference, because there is still a Condel score
        a.set_value('Conserved','')
        assert a.priority() == 3

        # Setting value to < 0.7 should make it not conserved any more, dropping it down to priority 2
        a.set_value('Condel','0.5')
        assert a.priority() == 2

        # Set a conservation value, but it should NOT change because there is a Condel score
        a.set_value('Conserved','437;Name=lod=80')
        assert a.priority() == 2


    def testNovelRare(self):

        # test: variant index 3 can only happen to a novel variant

        a = self.a
        a.set_value('ESP5400_ALL','0.001')
        a.set_value('1000g2010nov_ALL','0.001')
        a.set_value('Condel','0.9')

        assert a.priority() == 1
        
    def testRareTruncating(self):
        a = self.a
        a.ExonicFunc = "stopgain SNV"
        assert a.priority() == 4

        a.set_value('ESP5400_ALL','0.001')
        assert a.priority() == 1
        
    def testNovelTruncating(self):
        a = self.a
        a.ExonicFunc = "stopgain SNV"

        # Make it novel
        a.set_value('ESP5400_ALL','')
        a.set_value('1000g2010nov_ALL','')
        a.set_value('dbSNP138','')
        assert a.priority() == 4

    def testMissenseVaryRare(self):
        a = self.a
        a.set_value('ESP5400_ALL','0.0001')
        a.set_value('1000g2010nov_ALL','')
        a.set_value('Conserved','')
        a.set_value('Condel','')
        assert a.priority() == 2

 
if __name__ == '__main__':
    unittest.main()             
