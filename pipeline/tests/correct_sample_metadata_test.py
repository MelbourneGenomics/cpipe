
import unittest
import imp
import os
import random
import re
import sys
import io

sys.path.append('../scripts/')
import correct_sample_metadata_file

class CorrectMetadataTest(unittest.TestCase):

  def test_correct_column(self):
    before = '  G4: a1 , a2  G3: b1 , b2   '
    assert correct_sample_metadata_file.correct_column( before ) == 'G4:a1,a2 G3:b1,b2'      

  def test_empty(self):
   src = io.StringIO( 'h1\tPrioritised_Genes\th3\nx1\t\t' )
   target = io.StringIO()
   correct_sample_metadata_file.correct_metadata( src, target )
   assert target.getvalue() == 'h1\tPrioritised_Genes\th3\nx1\t\t\n'

  def test_correct_metadata(self):
   src = io.StringIO( 'h1\tPrioritised_Genes\th3\nv1\t G4: a1, a2\tv3' )
   target = io.StringIO()
   correct_sample_metadata_file.correct_metadata( src, target )
   assert target.getvalue() == 'h1\tPrioritised_Genes\th3\nv1\tG4:a1,a2\tv3\n'

