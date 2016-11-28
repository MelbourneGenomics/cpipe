import unittest
import imp
import os
import random
import re
import sys
import io

sys.path.append('../scripts/')
import update_pipeline_run_id

class UpdatePipelineRunIDTest(unittest.TestCase):

    def test_generate_id(self):
       fn = 'tmp%i' % random.randint(0, 1e9)
       with open( fn, 'w' ) as fh:
         fh.write( 'site' )
       
       new_id = update_pipeline_run_id.generate_new_id( fn )
       os.remove( fn ) 
       assert new_id == 'site_000000001'

    def test_update_samples(self):
       src = io.StringIO( 'h1\th2\nl1\tl2' )
       target = io.StringIO()
       new_id = update_pipeline_run_id.write( src, target, 'site_123456789' )
       assert target.getvalue() == 'Pipeline_Run_ID\th1\th2\nsite_123456789\tl1\tl2'
      

