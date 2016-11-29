import sys

import re
from typing.io import TextIO
from cpipe_utility import read_metadata
from manage_batch import FIELDS

IS_NUMERIC = {'dna concentration', 'dna quantity', 'dna quality', 'mean coverage'}
IS_ENUMERATION = {
    'sex': {'Male', 'Female', 'Unknown', 'other'},
    'sample type': {'Normal', 'Tumour'},
    'consanguinity': {'No', 'Yes', 'Suspected', 'Unknown'},
    'ethnicity': {'Unknown', 'European', 'African', 'Asian'}
}
IS_DATE = {'dna_date', 'capture_date', 'sequencing_date'}


def is_valid_numeric(field):
    """
        is numeric field valid
    """
    if len(field) > 0:
        try:
            float(field)
            return True
        except ValueError:
            return False
    return True


def is_valid_enumeration(field, allowed):
    """
        is enumerated field valid
    """
    return len(field) == 0 or field in allowed


def is_valid_date(field):
    """
        is date field valid
    """
    return len(field) == 0 or len(field) == 8 and field.isdigit()


def validate_metadata(metadata: TextIO):
    """
        Validate the input file according to the Melbourne Genomics metadata file format specification
    """
    metadata = read_metadata(metadata)
    warnings = []

    for (i, series) in metadata.iterrows():
        for (j, column) in enumerate(series):
            if re.search('^\w+', column):
                warnings.append('Sample {0} field {1} (column {2}) contains leading whitespace'.format(i, metadata.column, j))




'''
    headers = metadata.readline().strip().split('\t')
    idx = -1
    warnings = []
    field_filter = set()
    if fields is not None:
        for field in fields.split(','):
            field_filter.add(field.lower().strip())

    for idx, line in enumerate(sample_fh):
        fields = line.strip('\n').split('\t')
        for jdx, field in enumerate(fields):
            if len(field_filter) > 0 and headers[jdx].lower() not in field_filter:
                continue
            out.write("{0:>24}: {1}\n".format(headers[jdx], field))
            if field.startswith(' '):
                warnings.append(
                    'Sample {0} field {1} (column {2}) contains leading whitespace'.format(idx, headers[jdx], jdx))
            if field.endswith(' '):
                warnings.append(
                    'Sample {0} field {1} (column {2}) contains trailing whitespace'.format(idx, headers[jdx], jdx))
            if headers[jdx].lower() in IS_NUMERIC and not is_valid_numeric(field):
                warnings.append(
                    'Sample {0} field {1} (column {2}) cannot be "{3}": must be empty or a number'.format(idx,
                                                                                                          headers[jdx],
                                                                                                          jdx, field))
            if headers[jdx].lower() in IS_ENUMERATION and not is_valid_enumeration(field, IS_ENUMERATION[
                headers[jdx].lower()]):
                warnings.append(
                    'Sample {0} field {1} (column {2}) cannot be "{3}": must be empty or one of: {4}'.format(idx,
                                                                                                             headers[
                                                                                                                 jdx],
                                                                                                             jdx, field,
                                                                                                             ', '.join(
                                                                                                                 IS_ENUMERATION[
                                                                                                                     headers[
                                                                                                                         jdx].lower()])))
            if headers[jdx].lower() in IS_DATE and not is_valid_date(field):
                warnings.append(
                    'Sample {0} field {1} (column {2}) cannot be "{3}": must be empty or date (yyyymmdd)'.format(idx,
                                                                                                                 headers[
                                                                                                                     jdx],
                                                                                                                 jdx,
                                                                                                                 field))

    if idx == -1:
        err.write("ERROR: file only contains one line. Are you using Windows style line feeds?\n")
    for warning in warnings:
        err.write("WARNING: {0}\n".format(warning))
    if len(warnings) == 0:
        err.write("No warnings\n")
        '''
