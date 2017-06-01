#!/usr/bin/env python3
"""
Executable script designed to assist users in manipulating batches
"""
import sys

from cpipe.arg_validation import *
from cpipe.batch import Batch
import argparse


def setup_parser(parser: argparse.ArgumentParser):
    """
    Creates the argument parser for the manage batch script
    """

    subparsers = parser.add_subparsers(dest='command')
    parser.set_defaults(func=execute)
    parser.add_argument('--mgha', '-m', required=False, default=False,
                        help='Use MGHA-specific validation rules')

    # list command
    subparsers.add_parser('list', help='Lists the batches in the current Cpipe installation')

    # create command
    create_batch_parser = subparsers.add_parser('create',
                                                help='Creates a new batch, including data, metadata file and configuration file')
    create_batch_parser.add_argument('name', type=batch_name, help='The name for the new batch')
    create_batch_parser.add_argument('--data', '-d', required=True, help='The fastq files to add to the batch',
                                     nargs='+', type=path_with_ext(['.fastq', '.gz']))
    create_batch_parser.add_argument('--exome', '-e', required=True,
                                     help='A bed file indicating which regions are covered by the sequencing '
                                          'procedure', type=path_with_ext(['.bed']))
    create_batch_parser.add_argument('--profile', '-p', required=False,
                                     help='The analysis profile (gene list) to use for '
                                          'the analysis of this batch', type=profile, default='ALL')
    create_batch_parser.add_argument('--force', '-f', required=False, default=False, nargs='?', const=True,
                                     help='Replace an existing batch with'
                                          ' that name, if it already exists')
    create_batch_parser.add_argument('--mode', '-m', required=False, default='link',
                                     help='Either "copy", "link" or "move":'
                                          " the method used to put the data files into the batch directory")

    # edit command
    edit_parser = subparsers.add_parser('edit', help='Edit the metadata file for the chosen batch')
    edit_parser.add_argument('batch', type=existing_batch, help='The name of the batch whose metadata file you '
                                                                'want to edit')
    edit_parser.add_argument('--editor', '-e', type=editor, required=False,
                             help='The name of the executable you want to use to edit '
                                  'the metadata file using. Defaults to visidata, included with cpipe '
                                  '(https://github.com/saulpw/visidata)', default='vd')

    # view command
    view_parser = subparsers.add_parser('view',
                                        help='View the metadatafile for the chosen batch in a human-readable format')
    view_parser.add_argument('batch', type=existing_batch, help='The name of the batch whose metadata file you '
                                                                'want to view')
    # view_parser.add_argument('--sample', '-s', required=False,
    #                          help='The sample ID of the single sample you want to view the metadata for')

    # validate command
    validate_parser = subparsers.add_parser('check',
                                            help='Validate the metadata file for the chosen batch')
    validate_parser.add_argument('batch', type=existing_batch, help='The name of the batch whose metadata file you '
                                                                    'want to validate')

    # add_sample command
    add_sample_parser = subparsers.add_parser('add_sample',
                                              help='Add a sample to an existing metadata file')
    add_sample_parser.add_argument('batch', type=existing_batch, help='The name of the batch to which you'
                                                                      'want to add a sample')
    add_sample_parser.add_argument('data', nargs='+', help='The list of fastqs you want to add as samples to the batch')

    return parser


def execute(args: argparse.Namespace):
    if args.command == 'list':
        for batch in Batch.list_all():
            print(batch.name)
    elif args.command == 'create':
        Batch.create(args.name, args.data, args.exome, args.profile, force=args.force, mode=args.mode)
    elif args.command == 'edit':
        args.batch.metadata.edit()
    elif args.command == 'view':
        args.batch.metadata.view()
    elif args.command == 'check':
        warnings = args.batch.metadata.validate()
        if warnings:
            for warning in warnings:
                print(warning, file=sys.stderr)
            sys.exit(1)
        else:
            print(f'The metadata file for batch "{args.batch.name}" successfully passed the metadata check!')
    elif args.command == 'add_sample':
        args.batch.add_sample(args.data)
    else:
        raise ValueError('Unknown command')


def main():
    parser = argparse.ArgumentParser(description='Manage Cpipe batches and metadata files')
    setup_parser(parser)
    args = parser.parse_args()
    execute(args)


if __name__ == '__main__':
    main()
