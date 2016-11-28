
import unittest
import imp
import os
import random
import re
import sys
import io

sys.path.append('../scripts/')
import mark_batch_finished

class MarkBatchFinishedTest(unittest.TestCase):

  def test_simple(self):
    log = io.StringIO()
    mark_batch_finished.mark_batch_finished("/batch/dummy/analysis", log, read_only=False, move=False, dry=True)
    lines = log.getvalue().split('\n')
    assert len(lines) == 1 # finishes with blank line

  def test_full(self):
    log = io.StringIO()
    mark_batch_finished.mark_batch_finished("/batch/dummy/analysis", log, read_only=True, move=True, dry=True)
    lines = log.getvalue().split('\n')
    assert lines[0] == 'would execute: find "/batch/dummy" -path "*/.bpipe" -prune -o -name commandlog.txt -prune -o -exec chmod -w "{}" \\;'
    assert lines[1].startswith('would execute: mv /batch/dummy /batch/complete.dummy.')
    assert lines[2].startswith('would execute: ln -s /batch/complete.dummy.')
    assert len(lines) == 4 # finishes with blank line
