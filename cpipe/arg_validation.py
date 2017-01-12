import argparse
import shutil
from pathlib import Path

from cpipe.batch import Batch
from cpipe.design import Design


def existing_batch(batch_name):
    batch = Batch.find_by_name(batch_name)
    if batch:
        return batch
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
    design = Design.find_by_name(name)
    if design:
        return design
    else:
        raise argparse.ArgumentTypeError(
            'The profile must be an existing profile, namely a directory within the designs directory')


def editor(path):
    if shutil.which(path):
        return path
    else:
        raise argparse.ArgumentTypeError(
            'The editor you provide must be a valid executable accessible within your $PATH!')

__all__ = [existing_batch, batch_name, path_with_ext, profile, editor]
