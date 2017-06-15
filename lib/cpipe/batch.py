import os
import typing
from typing import List
import stat
from itertools import groupby
import subprocess
from pathlib import Path

from .paths import BATCHES
from .design import Design
from .metadata import Metadata
from . import pathlib_patches


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

    @classmethod
    def create(cls,
               batch_name: str,
               data: typing.Iterable[Path], 
               exome: str, 
               profile: Design = Design('ALL'),
               force: bool = False,
               mode: str = 'link',
               check_md5 = True,
               metadata: Path = None):
        """
        Creates a new batch and associated config files
        :param batch_name: The name for the new batch
        :param data: The path to the files to add to this batch
        :param exome: The exome target bed file
        :param profile: The design that specifies the regions to use in analysing this batch
        :param force: If a batch with this name already exists, delete it
        :param mode: 'link', 'copy', or 'move'
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

        # Add the metadata file if present
        if metadata is not None:
            cls.add_file(metadata, batch.metadata.path, mode='copy', check_md5=check_md5, force=True)

        # Write the config file
        with batch.config_file.open('w') as config_file:
            config_file.write(f'EXOME_TARGET="{exome}"')

        # Group fastqs into samples
        for id, fastqs in groupby(data, lambda x: x.name.split('_')[0]):

            fastqs = list(fastqs)

            # Move the data into the batch
            for fastq in fastqs:
                batch.add_fastq(fastq, mode=mode, force=force, check_md5=check_md5)

            # Update the metadata file if we don't already have one
            if metadata is None:
                batch.metadata.add_samples(fastqs, profile)

    def add_fastq(self, fastq: Path, mode: str = 'link', check_md5=True, force=False):
        """
        Physically adds a fastq file to the given batch using the method specified by mode
        :param fastq: The path to the fastq to add
        :param mode: Either 'link', 'copy', or 'move', indicating the method to be used to add the fastq to the batch
        """

        # Make the data subdir
        self.data.mkdir(exist_ok=True)

        # Move the data into the batch
        self.add_file(fastq, self.data, mode, force=force, check_md5=check_md5)

    @classmethod
    def add_file(cls, file: Path, dest: Path, mode: str = 'link', check_md5=True, force=False):
        """Add a file to a given directory using the given mode. Either 'link', 'move', or 'copy'"""

        # dest can be a filepath or a directory
        if dest.is_dir():
            target = dest / file.name
        else:
            target = dest

        # Make sure the destination file doesn't exist
        if target.exists():
            if force:
                target.unlink()
            else:
                raise Exception('The file "{target}" already exists. Check the --help for the command you ran to find a force option.')

        # Check any md5 sums that exist
        if check_md5:
            cls.check_md5(file)

        if mode == 'copy':
            file.copy(target)
        elif mode == 'link':
            file.symlink_from(target)
        elif mode == 'move':
            file.rename(target)

        # Fix permissions
        # target.chmod(stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)

    @staticmethod
    def check_md5(file: Path):
        suffix = file.suffix
        md5 = file.with_suffix(suffix + '.md5')
        if md5.exists():
            code = subprocess.run(['md5sum', '-c', md5], cwd=file.parent)
            try:
                code.check_returncode()
            except:
                raise Exception(f'md5sum of file {file} was unsuccessful. Check the integrity of the file. If you wish to turn this warning'
                        'off, check the --help for the command you just ran.')

    def add_samples(self, samples: List[Path], design: Design = Design('ALL'), mode: str = None):
        """
        Adds the samples to the batch directory and updates the metadata file
        :param samples:
        :param design:
        :param mode:
        :return:
        """

        for sample in samples:
            self.add_fastq(sample, mode)

        self.metadata.add_samples(samples, design)

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
        return Metadata(self.metadata_file, self)

    @property
    def metadata_file(self):
        return self.path / 'samples.txt'

    @property
    def config_file(self):
        return self.path / 'config.batch.groovy'

    @property
    def data(self):
        return self.path / 'data'

    @property
    def analysis(self):
        return self.path / 'analysis'

    def create_empty(self):
        """
        Creates all the important files and directories, but does not populate them with any files or content
        """
        self.path.mkdir(exist_ok=True)
        self.data.mkdir(exist_ok=True)
        self.analysis.mkdir(exist_ok=True)
        self.metadata.create_empty()
        self.config_file.touch(exist_ok=True)

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


__all__ = ["Batch"]
