import subprocess
import typing
from pathlib import Path
import pandas as pd

from .metadata_schema import get_schema
import cpipe
import cpipe.batch
import cpipe.design

class Metadata:
    def __init__(self, path: Path, batch: 'batch.Batch'):
        self.batch = batch
        self.path = path

    @property
    def data_frame(self):
        """
        Returns the metadata file as a pandas data frame
        """
        return cpipe.read_metadata(self.path)

    def create_empty(self):
        """
        Creates an empty metadata file, containing only the headers
        """
        schema = get_schema(is_mgha=False)
        names = '\t'.join([column.name for column in schema.columns])

        with self.path.open('w') as file:
            file.write(names)


    def view(self):  # sample: str = None):
        """
        Shows the sample metadata file using a spreadsheet viewer
        """
        subprocess.run(['vd', str(self.path), '--readonly'])

    def validate(self, is_mgha: bool = True):
        """
        Validate the input file according to the Melbourne Genomics metadata file format specification
        """

        metadata = self.data_frame
        schema = get_schema(is_mgha)
        return schema.validate(metadata)

    def edit(self, editor: str = 'vd', is_mgha: bool = False):
        """
        Edits the sample metadata file using a spreadsheet editor
        """
        subprocess.run([editor, str(self.path)])
        self.validate(is_mgha)

    def add_sample(self, sample: Path, design: 'design.Design' = None):
        """
        Adds a sample to the metadata file
        :param sample: The sample to add
        :param design: The name of the cohort to use for analysis. Defaults to the most common design currently in the
            metadata file
        """
        self.add_samples([sample], design)

    def add_samples(self, samples: typing.Iterable[Path], design: 'design.Design' = None):
        """
        Adds a number of samples to the metadata file
        :param samples: A list of samples to add
        :param design: The name of the cohort to use for analysis. Defaults to the most common design currently in the
            metadata file
        """

        # The sample ID is the text in the fastq filename before the first underscore
        ids = [sample.stem.split('_')[0] for sample in samples]
        id = ids[0]
        if ids.count(id) != len(ids):
            raise ValueError(
                'All fastqs from the same sample must have the same id (the text before the first underscore)')

        df = self.data_frame

        # Unless otherwise specified, the design is whatever is most common so far
        if not design:
            design = df['Cohort'].mode

        df.append({
            'Batch': self.batch.name,
            'Sample_ID': id,
            'Sex': 'Unknown',
            'Cohort': design.name,
            'Sample_Type': 'Normal',
            'Fastq_Files': ','.join([str(f.resolve()) for f in samples])
        }, ignore_index=True)\
            .to_csv(str(self.path), sep='\t', index=False)
