
import unittest
import imp
import os
import random
import re
import sys
import StringIO

sys.path.append('../scripts/')
import update_metadata

class UpdateMetadataTest(unittest.TestCase):

    def test_simple(self):
        log = StringIO.StringIO()
        sin = ("sample_id\tfoo", "1\t2")
        sout = StringIO.StringIO()
        update_metadata.update_metadata(sin, sout, log, "1", "foo", "bar")

        lines = sout.getvalue().split('\n')
        assert lines[0] == 'sample_id\tfoo1\tbar'
        assert len(lines) == 2 # finishes with blank line

    def test_no_sample(self):
        log = StringIO.StringIO()
        sin = ("sample_id\tfoo", "1\t2")
        sout = StringIO.StringIO()
        update_metadata.update_metadata(sin, sout, log, "2", "foo", "bar")

        lines = log.getvalue().split('\n')
        assert lines[0].startswith('ERROR')
        assert len(lines) == 2 # finishes with blank line
