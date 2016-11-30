#!/usr/bin/env python3
import typing

import pandas as pd
import re
import sys
import subprocess
from pathlib import Path
from typing import List
from itertools import groupby
import datetime

import pathlib_patches
from cpipe_utility import CONFIG_GROOVY_UTIL, CLASSPATH, BASE, BATCHES, DESIGNS, batch_dir, read_metadata
import cpipe_utility

DEFAULT_PRIORITY = '1'
FIELDS = ["Batch", "Sample_ID", "DNA_Tube_ID", "Sex", "DNA_Concentration", "DNA_Volume", "DNA_Quantity", "DNA_Quality",
          "DNA_Date", "Cohort", "Sample_Type", "Fastq_Files", "Prioritised_Genes", "Consanguinity", "Variants_File",
          "Pedigree_File", "Ethnicity", "VariantCall_Group", "Capture_Date", "Sequencing_Date", "Mean_Coverage",
          "Duplicate_Percentage", "Machine_ID", "DNA_Extraction_Lab", "Sequencing_Lab", "Exome_Capture",
          "Library_Preparation", "Barcode_Pool_Size", "Read_Type", "Machine_Type", "Sequencing_Chemistry",
          "Sequencing_Software", "Demultiplex_Software", "Hospital_Centre", "Sequencing_Contact", "Pipeline_Contact",
          "Notes", "Pipeline_Notes", "Analysis_Type"]
"""


def add_batch(batch_name, profile_name, exome_name, data_files, force, log):
    '''
        create a new batch
    '''
    batch_directory = batch_dir(batch_name)
    # check running from correct directory
    if not os.path.isdir(BATCHES):
        write_log(log, "ERROR: batches directory not found. Are you running this in the Cpipe root folder?")
        sys.exit(1)

    # make batch directory
    if not os.path.isdir(batch_directory):
        write_log(log, "Creating directory: {0}...".format(batch_directory))
        os.mkdir(batch_directory)
        write_log(log, "Creating directory: {0}: done".format(batch_directory))
    else:
        write_log(log, "Batch directory {0} exists".format(batch_directory))
        if os.path.exists(os.path.join(batch_directory, "samples.txt")):
            if force:
                write_log(log, "WARNING: samples.txt will be overwritten")
            else:
                write_log(log, "ERROR: samples.txt exists. Remove or specify --force to overwrite")
                sys.exit(1)

    # ensure profile exists
    profile_directory = os.path.join(DESIGNS, profile_name)
    profile = '{0}/{1}/{1}.genes.txt'.format(DESIGNS, profile_name)
    profile_bed = '{0}/{1}/{1}.bed'.format(DESIGNS, profile_name)
    if not os.path.isfile(profile):
        write_log(log, 'ERROR: profile "{0}" not found. Use manage_genelists.py to list or add new profiles'.format(
            profile_name))
        sys.exit(1)

    # configure exome target if specified
    # create_batch.sh behaviour:
    # - if no exome specified, look for target.bed
    # - look in designs for specified exome
    # - if not there, make one from the design gene list and write to initial-design.bed
    # - write exome_target to config.batch.groovy

    if exome_name is None:
        exome_name = '{0}.bed'.format(profile_name)

    # look in designs
    if os.path.exists(os.path.join(profile_directory, profile_name)):
        target_region = os.path.abspath(os.path.join(profile_directory, profile_name))
    # look at abspath
    elif os.path.exists(exome_name):
        target_region = os.path.abspath(exome_name)
    # generate the target regions
    else:
        target_region = os.path.join(profile_directory, "generated_target_regions.bed")
        write_log(log, 'Generating target region from design details and writing to {0}...'.format(target_region))
        run_command(
            "python {base}/pipeline/scripts/combine_target_regions.py --genefiles {0} --bedfiles {1} --exons {2} > {3}".format(
                profile, profile_bed, os.path.join(DESIGNS, "genelists/exons.bed"), target_region, base=BASE), log)
        write_log(log, 'Generating target region from design details and writing to {0}: done'.format(target_region))
        target_region = os.path.abspath(target_region)

    # write this to the batch
    target_file = os.path.join(batch_directory, "config.batch.groovy")
    write_log(log, 'writing {0}...'.format(target_file))
    with open(target_file, 'w') as target_fh:
        target_fh.write('EXOME_TARGET="{0}"\n'.format(target_region))
    write_log(log, 'writing {0}: done'.format(target_file))

    # generate samples.txt
    target_file = os.path.join(batch_directory, "samples.txt")
    if data_files is None:
        data_files = "data/*.fastq.gz"
    else:
        data_files = ' '.join(data_files)
    write_log(log, 'writing {0}...'.format(target_file))
    # for now we outsource sample file creation to the groovy
    run_command(
        '''
            set -e
            source {config_util}
            load_config
            cd {batches}/{batch_name}
            $GROOVY\
                -cp "{classpath}/*"\
                {base}/pipeline/scripts/files_to_sample_info.groovy\
                -batch {batch_name}\
                -disease {profile_name} {data_files}\
                > samples.txt
        '''.format(
            config_util=CONFIG_GROOVY_UTIL,
            batches=BATCHES, batch_name=batch_name, classpath=CLASSPATH,
            base=BASE,
            profile_name=profile_name,
            data_files=data_files
        ), log)
    write_log(log, 'writing {0}: done'.format(target_file))

    analysis_directory = os.path.join(batch_directory, "analysis")
    if not os.path.isdir(analysis_directory):
        write_log(log, "Creating directory: {0}...".format(analysis_directory))
        os.mkdir(analysis_directory)
        write_log(log, "Creating directory: {0}: done".format(analysis_directory))
    write_log(log,
              "To run an analysis: `./cpipe --batch {} run`".format(batch_name))


def add_sample(batch_name, profile_name, data_files, log):
    '''
        add to an existing sample
    '''
    batch_directory = batch_dir(batch_name)
    # generate samples.txt
    target_file = os.path.join(batch_directory, "samples.txt")
    if data_files is None:
        data_files = "data/*.fastq.gz"
    else:
        data_files = ' '.join(data_files)

    write_log(log, 'writing {0}...'.format(target_file))
    run_command(
        '''
            set -e
            source {config_util}
            load_config
            cd {batches}/{batch_name}
            $GROOVY\
                -cp "{classpath}/*"\
                {base}/pipeline/scripts/files_to_sample_info.groovy\
                -noheader\
                -disease {profile_name} {data_files}\
                >> samples.txt
        '''.format(
            config_util=CONFIG_GROOVY_UTIL,
            batch_name=batch_name,
            batches=BATCHES,
            classpath=CLASSPATH,
            base=BASE,
            profile_name=profile_name,
            data_files=data_files
        ), log)
    write_log(log, 'writing {0}: done'.format(target_file))


def show_batch(batch_name, out):
    '''
        show details of a single batch
    '''
    # batch exists?
    batch_directory = batch_dir(batch_name)
    if not os.path.isdir(batch_directory):
        out.write('Batch {0}: not found\n'.format(batch_name))
        sys.exit(0)
    # samples.txt
    samples = os.path.join(batch_directory, 'samples.txt')
    if os.path.exists(samples):
        lines = open(samples, 'r').readlines()
        out.write('BATCH {0} contains {1} sample(s)\n'.format(batch_name, len(lines) - 1))
        header = lines[0].strip('\n').split('\t')
        justify = max([len(x) for x in header])
        aggregate = collections.defaultdict(set)
        for line in lines[1:]:
            fields = line.strip('\n').split('\t')
            for idx, field in enumerate(fields):
                if len(field) > 0:
                    aggregate[header[idx]].add(field)
        # display aggregate of sample:
        for header in sorted(aggregate.keys()):
            out.write('{0}: {1}\n'.format(header.rjust(justify), ' '.join(aggregate[header])))
    else:
        out.write('samples.txt: not found\n'.format(batch_name))
        # could also show data, target regions

"""


def list_batches():
    df = cpipe_utility.list_batches()
    df.to_csv(sys.stdout, sep='\t', index=False)


def add_fastq(batch: Path, fastq: Path, mode: str = 'link'):
    """
    Adds a fastq file to the given batch using the method specified by mode
    :param batch:
    :param fastq:
    :param mode:
    :return:
    """

    # Make the data subdir
    data_dir = batch / 'data'
    data_dir.mkdir()

    # Move the data into the batch
    target_fastq = data_dir / Path(fastq).stem
    if mode == 'copy':
        fastq.copy(target_fastq)
    elif mode == 'link':
        fastq.symlink_from(target_fastq)
    elif mode == 'move':
        fastq.rename(target_fastq)


def add_sample_to_metadata(samples: List[Path], design: str, batch: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    # The sample ID is the text in the fastq filename before the first underscore
    ids = [sample.stem.split('_')[0] for sample in samples]
    id = ids[0]
    if ids.count(id) != len(ids):
        raise ValueError('All fastqs from the same sample must have the same id (the text before the first underscore)')

    return metadata.append({
        'Batch': batch.stem,
        'Sample_ID': id,
        'Sex': 'Unknown',
        'Cohort': design,
        'Sample_Type': 'Normal',
        'Fastq_Files': ','.join([str(f.resolve()) for f in samples])
    })


def create_batch(batch: Path, data: list, exome: str, profile: str, force: bool = False, mode: str = 'link'):
    # Handle an existing batch
    if batch.exists():
        if force:
            batch.rmtree()
        else:
            raise FileExistsError('Batch directory already exists! Use the --force flag to replace an existing batch')

    # Create the batch directory
    batch.mkdir()

    # Make the metadata file and open it
    metadata = batch / 'samples.txt'
    with metadata.open('w') as metadata_file:
        df = read_metadata(metadata_file)

        # Group fastqs into samples
        for (id, fastqs) in groupby(data, lambda x: x.split('_')[0]):

            # Move the data into the batch
            for fastq in fastqs:
                add_fastq(batch, fastq, mode=mode)

            # Update the metadata file
            add_sample_to_metadata(fastqs, profile, batch, df)

        # Write out the CSV
        df.to_csv(metadata_file, sep='\t')


def edit_batch(batch: Path, editor: str = 'editor'):
    metadata = batch / 'samples.txt'
    subprocess.run([editor, metadata])


def view_batch(batch: Path, sample: str = None):
    # Read the metadata file
    metadata = batch / 'samples.txt'
    df = read_metadata(metadata)

    # Subset the data frame if we only want one sample
    if sample:
        df = df[df['Sample_ID' == sample]]

    # Print out each row of the metadata file
    for (index, series) in df.iterrows():
        print(series)


def process_validation(series: pd.Series, message: str, error_list=None) -> (
        int, str, str):
    """
    Takes the given boolean series, and generates a warning for each item in the series that has a value of False
    :param series: A pandas boolean series. Each item should be True if it passed the validation, and False if it
    didn't.
    :param message: The message to print out, after the row and column numbers, indicating the validation error
    :param error_list: An optional array to add the validation errors to
    :return: Returns the array of warnings discovered using this validation
    """
    warnings = ['Row {} in column "{}" {}'.format(i + 1, series.name, message) for i, row in enumerate(series) if
                not row]

    if error_list is not None and error_list.extend:
        error_list.extend(warnings)

    return warnings


def validate_metadata(batch: Path):
    """
        Validate the input file according to the Melbourne Genomics metadata file format specification
    """
    metadata = read_metadata(batch / 'samples.txt', parse_num=False)
    warnings = []

    # Validate leading and trailing whitespace
    for column in metadata.columns:
        series = metadata[column]
        process_validation(~series.str.contains('^\s+'), 'has leading whitespace', warnings)
        process_validation(~series.str.contains('\s+$'), 'has trailing whitespace', warnings)

    # Validate Batch
    batches = cpipe_utility.list_batches()
    process_validation(~metadata['Batch'].isin(batches['Batch Name']), 'does not correspond to a valid batch directory',
                       warnings)

    # Validate Sample_ID
    process_validation(~metadata['Sample_ID'].str.contains('_'), 'contains an underscore', warnings)

    # Validate Sex
    sex_options = ['male', 'm', 'female', 'f', 'unknown']
    process_validation(metadata['Sex'].str.lower().isin(sex_options), 'must be one of {{{}}}'.format(', '.join(sex_options)),
                       warnings)

    # Validate numeric columns
    for column in ('DNA_Concentration', 'DNA_Volume', 'DNA_Quantity', 'DNA_Quality', 'Mean_Coverage', 'Duplicate_Percentage', 'Barcode_Pool_Size'):
        series = metadata[column]
        def valid_float(val):
            if len(val) == 0:
                return True
            try:
                float(val)
                return True
            except:
                return False

        process_validation(series.apply(valid_float), 'does not contain a valid numeric value', warnings)


    # Validate date columns
    for column in ('DNA_Date', 'Capture_Date', 'Sequencing_Date'):
        series = metadata[column]
        def valid_date(val):
            if len(val) == 0:
                return True
            try:
                datetime.datetime.strptime(val, '%Y%m%d')
                return True
            except:
                return False

        process_validation(series.apply(valid_date), 'does not contain a valid date in the format YYYYMMDD', warnings)



    '''
       Pedigree file column may be either a pedigree string, empty, import, individual or exclude (case insensitive)
    '''

    '''
        Gender may be Male, Female, Unknown or Other
    '''

    # TODO: Ensure all fastqs are in metadata file, and all metadata file fastqs exist

    '''
    for (i, series) in metadata.iterrows():
        for (j, column) in enumerate(series):
            if re.search('^\w+', column):
                warnings.append(
                    'Sample {0} field {1} (column {2}) contains leading whitespace'.format(i, metadata.column, j))
                    '''

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


'''

    if args.command == 'list':  # list all batches
        list_batches()
    else:
        if not args.batch:
            parser.print_help()
            sys.exit(1)
        if args.command == 'show_batch':  # show details of a batch
            show_batch(args.batch, out=sys.stdout)
        elif args.command == 'add_batch':  # add a new batch
            if not args.profile:
                write_log(sys.stderr, "No profile specified; defaulting to ALL (entire human exome)")
                add_batch(args.batch, 'ALL', args.exome, args.data, args.force, log=sys.stderr)
            else:
                add_batch(args.batch, args.profile, args.exome, args.data, args.force, log=sys.stderr)
        elif args.command == 'add_sample':  # add additional samples to batch
            add_sample(args.batch, args.profile, args.data, log=sys.stderr)

            '''
