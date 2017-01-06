import os
import typing
import subprocess
from pathlib import Path
from typing import List
from itertools import groupby
import pandas as pd
import sys

from cpipe_util.paths import CONFIG_GROOVY_UTIL, CLASSPATH, BASE, BATCHES, DESIGNS
from cpipe_util import read_metadata, Batch, Design, pathlib_patches
import cpipe_util
from manage_batch.schema import get_schema

from pathlib import Path
from .paths import BATCHES
from cpipe_util import read_metadata, pathlib_patches, Design

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
    def create(batch_name: str, data: typing.Iterable[Path], exome: str, profile: Design, force: bool = False, mode: str = 'link'):
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
                raise FileExistsError('Batch directory already exists! Use the --force flag to replace an existing batch')

        # Create the batch directory
        batch.create()

        # Write the config file
        with batch.config.open('w') as config_file:
            config_file.write(f'EXOME_TARGET={exome}')

        # Make the metadata file and open it
        with batch.metadata.open('w') as metadata_file:
            df = pd.DataFrame()

            # Group fastqs into samples
            for id, fastqs in groupby(data, lambda x: x.split('_')[0]):

                # Move the data into the batch
                for fastq in fastqs:
                    batch.add_fastq(batch, fastq, mode=mode)

                # Update the metadata file
                batch.add_sample_to_metadata(fastqs, profile, batch, df)

            # Write out the CSV
            df.to_csv(metadata_file, sep='\t')

    def add_fastq(self, fastq: Path, mode: str = 'link'):
        """
        Adds a fastq file to the given batch using the method specified by mode
        :param batch:
        :param fastq:
        :param mode:
        :return:
        """

        # Make the data subdir
        data_dir = self.path / 'data'
        data_dir.mkdir()

        # Move the data into the batch
        target_fastq = data_dir / Path(fastq).stem
        if mode == 'copy':
            fastq.copy(target_fastq)
        elif mode == 'link':
            fastq.symlink_from(target_fastq)
        elif mode == 'move':
            fastq.rename(target_fastq)

    def add_sample_to_metadata(self, samples: List[Path], design: Design, metadata: pd.DataFrame) -> pd.DataFrame:

        # The sample ID is the text in the fastq filename before the first underscore
        ids = [sample.stem.split('_')[0] for sample in samples]
        id = ids[0]
        if ids.count(id) != len(ids):
            raise ValueError('All fastqs from the same sample must have the same id (the text before the first underscore)')

        return metadata.append({
            'Batch': self.name,
            'Sample_ID': id,
            'Sex': 'Unknown',
            'Cohort': design.name,
            'Sample_Type': 'Normal',
            'Fastq_Files': ','.join([str(f.resolve()) for f in samples])
        })



    def edit_batch(self, editor: str = 'editor', is_mgha: bool=False):
        subprocess.run([editor, str(self.metadata)])
        self.validate_metadata(is_mgha)


    def view_batch(self): #sample: str = None):
        subprocess.run(['vd', str(self.metadata)])

        # # Read the metadata file
        # metadata = batch.path / 'samples.txt'
        # df = cpipe_util.read_metadata(metadata)
        #
        # # Subset the data frame if we only want one sample
        # if sample:
        #     df = df[df['Sample_ID' == sample]]
        #
        # # Print out each row of the metadata file
        # for (index, series) in df.iterrows():
        #     print(series)


    def validate_metadata(self, is_mgha: bool = True):
        """
            Validate the input file according to the Melbourne Genomics metadata file format specification
        """

        metadata = read_metadata(self.path / 'samples.txt', parse_num=False)
        schema = get_schema(is_mgha)
        warnings = schema.validate(metadata)

        if warnings:
            for warning in warnings:
                print(warning, file=sys.stderr)
            sys.exit(1)
        else:
            print(f'The metadata file for batch "{self}" successfully passed the metadata check!')


    def add_sample(self, samples: List[Path]):
        metadata = self.metadata_df

        # The design is whatever is most common so far
        design = metadata['Cohort'].mode

        # Update the metadata file and save it
        self.add_sample_to_metadata(samples, design, metadata)
        metadata.to_csv(self.metadata)

    def delete(self):
        self.path.rmtree()

    @property
    def name(self):
        return str(self.path.relative_to(BATCHES))

    @property
    def metadata(self):
        return self.path / 'samples.txt'

    @property
    def metadata_df(self):
        """
        Returns a pandas DataFrame of the metadata file
        :return:
        """
        return read_metadata(self.metadata)

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
        batches = []

        for root, dirs, files in os.walk(str(BATCHES)):
            if 'samples.txt' in files:
                batches.append(Batch(root))

        # Sort by batch name
        return sorted(batches, key=lambda x: x.name)
