import unittest
from cpipe_util import Batch, paths, Design

class Batch_Test(unittest.TestCase):
    batch = Batch(name='test')
    fastqs = list((paths.BASE / 'pipeline/tests/detect_mutations_test/data').iterdir())
    source_fastq = fastqs[0]
    dest_fastq = batch.data / source_fastq.name


    def create_test_batch(self):
        Batch.create(
            self.batch.name,
            self.fastqs,
            paths.DESIGNS / 'genelists/exons.bed',
            Design('ALL')
        )

    def delete_test_batch(self):
        self.batch.delete()

    def test_create_delete_batch(self):
        """
        Creates a new batch called 'test', checks that all the relevant files are created, deletes it, and check that
        all the files were deleted
        """
        self.create_test_batch()

        # Assert
        self.assertTrue(self.batch.path.exists())
        self.assertTrue(self.batch.config_file.exists())
        self.assertTrue(self.batch.metadata.path.exists())
        self.assertTrue(len(self.batch.metadata.data_frame) > 0)
        self.assertTrue(self.batch.config_file.stat().st_size > 0)
        self.assertTrue(self.batch.metadata.path.stat().st_size > 0)

        # Cleanup data
        self.delete_test_batch()

        self.assertFalse(self.batch.path.exists())
        self.assertFalse(self.batch.config_file.exists())
        self.assertFalse(self.batch.metadata.path.exists())


    def test_add_sample_copy(self):
        """
        Checks that add_sample works, and copies correctly
        """
        self.batch.create_empty()
        self.batch.add_fastq(self.source_fastq, mode='copy')

        # Asserts
        self.assertTrue(self.source_fastq.exists())
        self.assertFalse(self.source_fastq.is_symlink())
        self.assertTrue(self.dest_fastq.exists())
        self.assertFalse(self.dest_fastq.is_symlink())

        self.delete_test_batch()

    def test_add_sample_link(self):
        """
        Checks that add_sample works, and links correctly
        """

        self.batch.create_empty()
        self.batch.add_fastq(self.source_fastq, mode='link')

        # Asserts
        self.assertTrue(self.source_fastq.exists())
        self.assertFalse(self.source_fastq.is_symlink())
        self.assertTrue(self.dest_fastq.exists())
        self.assertTrue(self.dest_fastq.is_symlink())

        self.delete_test_batch()

    def test_add_sample_move(self):
        """
        Checks that add_sample works, and moves correctly
        """
        self.batch.create_empty()
        self.batch.add_fastq(self.source_fastq, mode='move')

        # Asserts
        self.assertFalse(self.source_fastq.exists())
        self.assertTrue(self.dest_fastq.exists())
        self.assertFalse(self.dest_fastq.is_symlink())

        # Move the file back
        self.dest_fastq.rename(self.source_fastq)
        self.delete_test_batch()
