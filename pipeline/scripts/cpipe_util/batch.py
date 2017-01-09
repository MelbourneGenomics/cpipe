import os
import typing
from typing import List
from itertools import groupby
import pandas as pd

from pathlib import Path
from .paths import BATCHES
from .design import Design
from .metadata import Metadata


class Batch:

    def __init__(self, directory=None, name=None):
        if directory and name:
            raise ValueError("You can't specify a directory and a name for the batch. Use one or either")
        elif directory:
            self.path = Path(directory).resolve()
        elif name:
            self.path = BATCHES / name
        else:
            raise ValueError("You must specify a directory or a name for the batch. Use one or either")

    @staticmethod
    def create(batch_name: str, data: typing.Iterable[Path], exome: str, profile: Design, force: bool = False,
               mode: str = 'link'):
        """
        Creates a new batch and associated config files
        :param batch_name: The name for the new batch
        :param data: The path to the files to add to this batch
        :param exome: The exome target bed file
        :param profile: The design that specifies the regions to use in analysing this batch
        :param force: If a batch with this name already exists, delete it
        :param mode: 'link', 'copy', or ''
        :return:
        """
        batch = Batch(BATCHES / batch_name)

        # Handle an existing batch
        if batch.path.exists():
            if force:
                batch.path.rmtree()
            else:
                raise FileExistsError(
                    'Batch directory already exists! Use the --force flag to replace an existing batch')

        # Create the batch directory
        batch.create_empty()

        # Write the config file
        with batch.config.open('w') as config_file:
            config_file.write(f'EXOME_TARGET={exome}')

        # Make the metadata file and open it
        with batch.metadata.path.open('w') as metadata_file:
            df = pd.DataFrame()

            # Group fastqs into samples
            for id, fastqs in groupby(data, lambda x: x.split('_')[0]):

                # Move the data into the batch
                for fastq in fastqs:
                    batch.add_fastq(batch, fastq, mode=mode)

                # Update the metadata file
                batch.metadata.add_sample(fastqs, profile, batch, df)

            # Write out the CSV
            df.to_csv(metadata_file, sep='\t')

    def add_fastq(self, fastq: Path, mode: str = 'link'):
        """
        Physically adds a fastq file to the given batch using the method specified by mode
        :param fastq: The path to the fastq to add
        :param mode: Either 'link', 'copy', or 'move', indicating the method to be used to add the fastq to the batch
        """

        # Make the data subdir
        self.data.mkdir()

        # Move the data into the batch
        target_fastq = self.data / Path(fastq).stem
        if mode == 'copy':
            fastq.copy(target_fastq)
        elif mode == 'link':
            fastq.symlink_from(target_fastq)
        elif mode == 'move':
            fastq.rename(target_fastq)


    def add_samples(self, samples: List[Path], design: Design, mode: str = None):
        """
        Adds the samples to the batch directory and updates the metadata file
        :param samples:
        :param design:
        :param mode:
        :return:
        """

        for sample in samples:
            self.add_fastq(sample, mode)

        self.add_samples(samples, design)

    def delete(self):
        """
        Deletes the entire batch directory and all its contents
        :return:
        """
        self.path.rmtree()

    @property
    def name(self):
        return str(self.path.relative_to(BATCHES))

    @property
    def metadata(self) -> Metadata:
        return Metadata(self.path / 'samples.txt', self)

    @property
    def config(self):
        return self.path / 'config.batch.groovy'

    @property
    def data(self):
        return self.path / 'data'

    def create_empty(self):
        """
        Creates all the important files and directories, but does not populate them with any files or content
        """
        self.path.mkdir(exist_ok=True)
        self.data.mkdir(exist_ok=True)
        self.metadata.touch(exist_ok=True)
        self.config.touch(exist_ok=True)

    @classmethod
    def find_by_name(cls, name: str) -> 'Batch':
        """
        Returns a Batch with the given name, or None if it does not exist
        :param name: A batch name
        """
        found = [batch for batch in cls.list_all() if batch.name == name]
        return found[0] if found else None

    @staticmethod
    def list_all() -> typing.Iterable['Batch']:
        """
            Prints the name of all batches that contain a samples.txt file
        """

        # Find all directories that contain a samples.txt and add them to a list
        batches = [Batch(root) for (root, dirs, files) in os.walk(str(BATCHES)) if 'samples.txt' in files]

        # Sort by batch name
        return sorted(batches, key=lambda x: x.name)

    def __str__(self):
        return self.name

    def __fspath__(self):
        return self.path
