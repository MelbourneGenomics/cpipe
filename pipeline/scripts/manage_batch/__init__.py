#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path
from typing import List
from itertools import groupby
import pandas as pd

from cpipe_util.paths import CONFIG_GROOVY_UTIL, CLASSPATH, BASE, BATCHES, DESIGNS
from cpipe_util import read_metadata, Batch, Design, pathlib_patches
import cpipe_util
from manage_batch.schema import get_schema

def list_batches():
    df = cpipe_util.list_batches()
    df.to_csv(sys.stdout, sep='\t', index=False)


def add_fastq(batch: Batch, fastq: Path, mode: str = 'link'):
    """
    Adds a fastq file to the given batch using the method specified by mode
    :param batch:
    :param fastq:
    :param mode:
    :return:
    """

    # Make the data subdir
    data_dir = batch.path / 'data'
    data_dir.mkdir()

    # Move the data into the batch
    target_fastq = data_dir / Path(fastq).stem
    if mode == 'copy':
        fastq.copy(target_fastq)
    elif mode == 'link':
        fastq.symlink_from(target_fastq)
    elif mode == 'move':
        fastq.rename(target_fastq)


def add_sample_to_metadata(samples: List[Path], design: Design, batch: Batch, metadata: pd.DataFrame) -> pd.DataFrame:

    # The sample ID is the text in the fastq filename before the first underscore
    ids = [sample.stem.split('_')[0] for sample in samples]
    id = ids[0]
    if ids.count(id) != len(ids):
        raise ValueError('All fastqs from the same sample must have the same id (the text before the first underscore)')

    return metadata.append({
        'Batch': batch.name,
        'Sample_ID': id,
        'Sex': 'Unknown',
        'Cohort': design.name,
        'Sample_Type': 'Normal',
        'Fastq_Files': ','.join([str(f.resolve()) for f in samples])
    })


def create_batch(batch_name: str, data: list, exome: str, profile: str, force: bool = False, mode: str = 'link'):
    """
    Creates a new batch and associated config files
    :param batch:
    :param data:
    :param exome:
    :param profile:
    :param force:
    :param mode:
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
    batch.path.mkdir()

    # Make the metadata file and open it
    metadata = batch.path / 'samples.txt'
    with metadata.open('w') as metadata_file:
        df = cpipe_util.read_metadata(metadata_file)

        # Group fastqs into samples
        for id, fastqs in groupby(data, lambda x: x.split('_')[0]):

            # Move the data into the batch
            for fastq in fastqs:
                add_fastq(batch, fastq, mode=mode)

            # Update the metadata file
            add_sample_to_metadata(fastqs, profile, batch, df)

        # Write out the CSV
        df.to_csv(metadata_file, sep='\t')


def edit_batch(batch: Batch, editor: str = 'editor', is_mgha: bool=False):
    metadata = batch.path / 'samples.txt'
    subprocess.run([editor, str(metadata)])
    validate_metadata(batch, is_mgha)


def view_batch(batch: Batch, sample: str = None):
    # Read the metadata file
    metadata = batch.path / 'samples.txt'
    df = cpipe_util.read_metadata(metadata)

    # Subset the data frame if we only want one sample
    if sample:
        df = df[df['Sample_ID' == sample]]

    # Print out each row of the metadata file
    for (index, series) in df.iterrows():
        print(series)


def validate_metadata(batch: Batch, is_mgha: bool = True):
    """
        Validate the input file according to the Melbourne Genomics metadata file format specification
    """

    metadata = read_metadata(batch.path / 'samples.txt', parse_num=False)
    schema = get_schema(is_mgha)
    warnings = schema.validate(metadata)

    if warnings:
        for warning in warnings:
            print(warning, file=sys.stderr)
        sys.exit(1)
    else:
        print('The metadata file for batch "{}" successfully passed the metadata check!'.format(batch))


def add_sample(batch: Path, samples: List[Path]):
    # Find the metadata file and read it
    metadata_file = batch / 'samples.txt'
    metadata = read_metadata(metadata_file)

    # The design is whatever is most common so far
    design = metadata['Cohort'].mode

    # Update the metadata file and save it
    add_sample_to_metadata(samples, design, batch, metadata)
    metadata.to_csv(metadata_file)

