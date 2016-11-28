import os
import unittest
import sys
import traceback
import io
import shutil
import argparse
from contextlib import contextmanager

sys.path.append('../scripts/')

from manage_batch import add_batch, show_batches, show_batch, add_sample, create_parser
from cpipe_utility import BASE


@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr

class ManageBatchTest(unittest.TestCase):
    test_batch_name = 'test'
    test_data_dir = os.path.join(BASE, 'pipeline/tests/detect_mutations_test/data')
    test_data_files = [os.path.join(test_data_dir, f) for f in os.listdir(test_data_dir)]
    test_batch_dir = os.path.join(BASE, 'batches', test_batch_name)
    batch_groovy = os.path.join(test_batch_dir, 'config.batch.groovy')
    sample_metadata = os.path.join(test_batch_dir, 'samples.txt')

    def create_test_batch(self, batch_name=None):
        batch_name = batch_name or self.test_batch_name
        with open(os.devnull, "w") as devnull:
            add_batch(batch_name, 'ALL', None, self.test_data_files, None, log=devnull)

    def delete_test_batch(self, batch_name=None):
        batch_name = batch_name or self.test_batch_name
        batch_dir = os.path.join(BASE, 'batches', batch_name)
        shutil.rmtree(batch_dir)

    def test_add_batch(self):
        """
        Creates a new batch called 'test'
        :return:
        """
        try:
            # Run the command itself
            with suppress_stdout():
                self.create_test_batch()

            # Assert
            assert (os.path.exists(self.batch_groovy))
            assert (os.path.exists(self.sample_metadata))
            assert (os.stat(self.batch_groovy).st_size > 0)
            assert (os.stat(self.sample_metadata).st_size > 0)

            # Cleanup data
            self.delete_test_batch()

        except Exception:
            self.fail('Failed with error: ' + traceback.format_exc())

    def test_show_batches(self):
        stream = io.StringIO()
        with suppress_stdout():
            show_batches(stream)
        assert (len(stream.getvalue()) > 0)

    def test_show_batch(self):
        self.create_test_batch()
        stream = io.StringIO()
        with suppress_stdout():
            show_batch(self.test_batch_name, stream)
        assert (len(stream.getvalue()) > 0)
        self.delete_test_batch()

    def test_create_parser(self):
        parser = create_parser()

        with self.assertRaises(SystemExit):
            parser.parse_args(['add_batch', '--profile', 'ALL', '--batch', 'batch'])
        with self.assertRaises(SystemExit):
            parser.parse_args(['add_batch', '--profile', 'ALL', '--batch', 'sample_id'])

    def test_add_sample(self):
        self.create_test_batch()
        initial_size = os.stat(self.sample_metadata).st_size
        with suppress_stdout(), open(os.devnull, "w") as devnull:
            add_sample(self.test_batch_name, 'ALL', self.test_data_files, log=devnull)
        final_size = os.stat(self.sample_metadata).st_size
        self.delete_test_batch()
        assert(final_size > initial_size)
