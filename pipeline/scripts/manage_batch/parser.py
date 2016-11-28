import argparse
import shutil
from pathlib import Path
from ..cpipe_utility import DESIGNS, DESIGN_LIST, list_batches


def existing_batch(val):
    if val in list_batches()['Batch Name']:
        return val
    else:
        raise argparse.ArgumentTypeError('The batch must be an existing batch with a metadata file')


def batch_name(val):
    if not isinstance(val, str):
        raise argparse.ArgumentTypeError('The batch name must be provided as a string')
    batch = str(val)
    if batch.lower().startswith(('batch', 'sample_id')):
        raise argparse.ArgumentTypeError('The batch name cannot start with "batch" or "sample_id" due to a bug in '
                                         'the Bpipe SampleInfo parser.')
    return batch


def path_with_ext(exts):
    exts = set(exts)

    def func(path):
        path = Path(path)
        if set(path.suffixes).issuperset(exts) and path.exists():
            return path
        else:
            raise argparse.ArgumentTypeError('The file must exist, and have a {} file extension'.format(exts))

    return func


def profile(name):
    if name in DESIGN_LIST:
        return name
    else:
        raise argparse.ArgumentTypeError(
            'The profile must be an existing profile, namely a directory within the designs directory')


def editor(path):
    if shutil.which(path):
        return path
    else:
        raise argparse.ArgumentTypeError(
            'The editor you provide must be a valid executable accessible within your $PATH!')


def create_parser():
    parser = argparse.ArgumentParser(description='Manage Cpipe batches and metadata files')
    subparsers = parser.add_subparsers()

    # list command
    list_parser = subparsers.add_parser('list', help='Lists the batches in the current Cpipe installation')

    # create command
    create_batch_parser = subparsers.add_parser('create',
                                                help='Creates a new batch, including data, metadata file and configuration file')
    create_batch_parser.add_argument('name', type=batch_name, required=True, help='The name for the new batch')
    create_batch_parser.add_argument('--data', '-d', required=True, help='The fastq files to add to the batch',
                                     nargs='+', type=path_with_ext(['fastq', '.gz']))
    create_batch_parser.add_argument('--exome', '-e', required=True,
                                     help='A bed file indicating which regions are covered by the sequencing '
                                          'procedure', type=path_with_ext(['bed']))
    create_batch_parser.add_argument('--profile', '-p', required=False,
                                     help='The analysis profile (gene list) to use for '
                                          'the analysis of this batch', type=profile, default='ALL')
    create_batch_parser.add_argument('--force', '-f', required=False, default=False, help='Replace an existing batch with'
                                     ' that name, if it already exists')

    # edit command
    edit_parser = subparsers.add_parser('edit', help='Edit the metadata file for the chosen batch')
    edit_parser.add_argument('batch', type=existing_batch, required=True,
                             help='The name of the batch whose metadata file you '
                                  'want to edit')
    edit_parser.add_argument('--editor', '-e', type=editor, required=False,
                             help='The name of the executable you want to use to edit '
                                  'the metadata file using', default='editor')

    # view command
    view_parser = subparsers.add_parser('view',
                                        help='View the metadatafile for the chosen batch in a human-readable format')
    view_parser.add_argument('batch', type=existing_batch, required=True,
                             help='The name of the batch whose metadata file you '
                                  'want to view')
    view_parser.add_argument('--sample', '-s', required=False,
                             help='The sample ID of the single sample you want to view the metadata for')

    # validate command
    validate_parser = subparsers.add_parser('validate', help='Validate the metadata file for the chosen batch')
    validate_parser.add_argument('batch', type=existing_batch, required=True,
                                 help='The name of the batch whose metadata file you '
                                      'want to view')

    # add_sample command
    add_sample_parser = subparsers.add_parser('add_sample', help='Add a sample to an existing metadata file')
    add_sample_parser.add_argument('batch', type=existing_batch, required=True,
                                   help='The name of the batch to which you'
                                        'want to add a sample')
    add_sample_parser.add_argument('samples', required=True, nargs='+',
                                   help='The list of fastqs you want to add as samples to the batch')

    return parser
