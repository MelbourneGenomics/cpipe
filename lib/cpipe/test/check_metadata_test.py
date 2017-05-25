
import unittest
import imp
import os
import random
import re
import sys
import io

from cpipe.scripts import check_metadata

class CheckMetadataTest(unittest.TestCase):

  def test_simple(self):
    src = io.StringIO('Batch\nNA12878')
    out = io.StringIO()
    err = io.StringIO()
    check_metadata.validate(src, out, err)
    lines = err.getvalue().split('\n')
    assert lines[0] == "No warnings"
    assert len(lines) == 2 # finishes with blank line


  def test_validate_empty(self):
    src = io.StringIO('')
    out = io.StringIO()
    err = io.StringIO()
    check_metadata.validate(src, out, err)
    lines = err.getvalue().split('\n')
    assert lines[0] == "ERROR: file only contains one line. Are you using Windows style line feeds?"

