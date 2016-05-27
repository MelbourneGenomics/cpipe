
import unittest
import imp
import os
import random
import re
import sys
import StringIO

sys.path.append('../scripts/')
import check_metadata

class CheckMetadataTest(unittest.TestCase):

  def test_validate_empty(self):
    src = StringIO.StringIO('')
    out = StringIO.StringIO()
    err = StringIO.StringIO()
    check_metadata.validate(src, out, err)
    lines = err.getvalue().split('\n')
    assert lines[0] == "ERROR: file only contains one line. Are you using Windows style line feeds?"

