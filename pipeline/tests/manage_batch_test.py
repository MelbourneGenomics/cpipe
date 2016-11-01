import os
import unittest
import sys

from manage_batch import add_batch
from cpipe_utilities import BASE

sys.path.append('../scripts/')

class ManageBatchTest(unittest.TestCase):
    def test_add_batch(self):
        """
        Creates a new batch called 'test'
        :return:
        """
        try:
            add_batch('test', None, None, os.listdir(os.path.join(BASE, 'batches/mutation_detection/data')), False, sys.stdout)
        except Exception as x:
            self.fail('Failed with error: ' + x.strerror)
