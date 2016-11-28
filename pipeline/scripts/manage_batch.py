#!/usr/bin/env python3

import argparse
import collections
import datetime
import glob
import pandas as pd
import os
import os.path
import fnmatch
import sys
import subprocess

from cpipe_utility import CONFIG_GROOVY_UTIL, CLASSPATH, BASE, BATCHES, DESIGNS, batch_dir

DEFAULT_PRIORITY = '1'
FIELDS = ["Batch", "Sample_ID", "DNA_Tube_ID", "Sex", "DNA_Concentration", "DNA_Volume", "DNA_Quantity", "DNA_Quality",
          "DNA_Date", "Cohort", "Sample_Type", "Fastq_Files", "Prioritised_Genes", "Consanguinity", "Variants_File",
          "Pedigree_File", "Ethnicity", "VariantCall_Group", "Capture_Date", "Sequencing_Date", "Mean_Coverage",
          "Duplicate_Percentage", "Machine_ID", "DNA_Extraction_Lab", "Sequencing_Lab", "Exome_Capture",
          "Library_Preparation", "Barcode_Pool_Size", "Read_Type", "Machine_Type", "Sequencing_Chemistry",
          "Sequencing_Software", "Demultiplex_Software", "Hospital_Centre", "Sequencing_Contact", "Pipeline_Contact",
          "Notes", "Pipeline_Notes", "Analysis_Type"]


def list_batches(out):
    '''
        Prints the name of all batches that contain a samples.txt file
    '''

    # Find all directories that contain a samples.txt and add them to a list
    df = pd.DataFrame(columns=('Batch Name', 'Batch Path'))
    for root, dirs, files in os.walk(BATCHES):
        if 'samples.txt' in files:
            batch_name = os.path.basename(root)
            full_path = os.path.abspath(root)
            df = df.append({'Batch Name': batch_name, 'Batch Path': full_path}, ignore_index=True)

    # Sort them alphabetically by their batch name
    df = df.sort_values(by='Batch Name')

    # Print to stdout
    df.to_csv(out, sep='\t', index=False)

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
            batches=BATCHES,
            batch_name=batch_name,
            classpath=CLASSPATH,
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


def batch_name(val):
    if not isinstance(val, str):
        raise argparse.ArgumentTypeError('The batch name must be provided as a string')
    batch = str(val)
    if batch.lower().startswith(('batch', 'sample_id')):
        raise argparse.ArgumentTypeError('The batch name cannot start with "batch" or "sample_id" due to a bug in '
                                         'the Bpipe SampleInfo parser.')
    return batch

def fastq_path(path):
    (base, ext) = os.path.splitext(path)
    if os.path.exists(path) and ext == 'gz' and base.endswith('fastq'):
        return path
    else:
        raise argparse.ArgumentTypeError('The fastq file must exist, and end in a ".fastq.gz" file extension')

def bed_path(path):
    (base, ext) = os.path.splitext(path)
    if os.path.exists(path) and ext == 'bed':
        return path
    else:
        raise argparse.ArgumentTypeError('The bed file must exist, and end in a ".bed" file extension')


def profile(name):
     if


def create_parser():
    parser = argparse.ArgumentParser(description='Manage Cpipe batches and metadata files')
    subparsers = parser.add_subparsers()

    # list command
    list_parser = subparsers.add_parser('list', help='Lists the batches in the current Cpipe installation')

    # create command
    create_batch_parser = subparsers.add_parser('create', help='Creates a new batch, including data, metadata file and configuration file')
    create_batch_parser.add_argument('name', type=batch_name, required=True, help='The name for the new batch')
    create_batch_parser.add_argument('--data', '-d', required=True, help='The fastq files to add to the batch', nargs='+', type=fastq_path)
    create_batch_parser.add_argument('--exome', '-e', required=True, help='A bed file indicating which regions are covered by the sequencing '
                                                                          'procedure', type=bed_path)
    create_batch_parser.add_argument('--profile', '-p', required=True, help='The analysis profile (genelist) to use for '
                                                                            'the analysis of this batch', type=fastq_path)

    parser.add_argument('command', help='command to execute',
                        choices=['list', 'create', 'edit_metadata', 'view_metadata', 'validate_metadata', 'add_sample'])
    parser.add_argument('--batch', required=False, help='batch name', type=validate_batch)
    parser.add_argument('--profile', required=False, help='analysis profile')
    parser.add_argument('--exome', required=False, help='target regions')
    parser.add_argument('--data', required=False, nargs='*', help='fastq files (relative to batch directory)')
    parser.add_argument('--force', action="store_true", required=False, help='force action')
    return parser


def main():
    '''
        parse command line
    '''

    parser = create_parser()
    args = parser.parse_args()

    if args.command == 'list':  # list all batches
        list_batches(out=sys.stdout)
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


if __name__ == '__main__':
    main()
