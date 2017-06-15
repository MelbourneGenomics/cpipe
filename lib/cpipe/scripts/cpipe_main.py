#!/usr/bin/env python3
import argparse
import unittest
from subprocess import check_output
from os import environ

import subprocess
from cpipe import arg_validation
from collections.abc import Sequence
from doit.doit_cmd import DoitMain
from cpipe.scripts import manage_batch, manage_genelists, check_java as cj
from cpipe import paths


def run(cmd, **kwargs):
    if isinstance(cmd, str):
        kwargs['shell'] = True
    if 'cwd' not in kwargs:
        kwargs['cwd'] = environ.get('CPIPE_ROOT')

    subprocess.run(cmd, **kwargs)


def check_java(args):
    """Run the java check if necessary"""
    if not args.no_java_check:
        cj.check_java()


def run_pipeline(args):
    check_java(args)
    args.batch.analysis.mkdir(exist_ok=True, parents=True)
    run(
        f'bpipe run {" ".join(args.bpipe_opts)} {paths.BASE}/pipeline/pipeline.groovy {args.batch.metadata_file}',
        cwd=args.batch.analysis
    )


def test_pipeline(args):
    check_java(args)
    testsuite = unittest.TestLoader().discover('cpipe.test', pattern='*')
    unittest.TextTestRunner(verbosity=1).run(testsuite)

    run('run_tests detect_mutations_test')


def stop_pipeline(args):
    check_java(args)
    run(
        'bpipe stop',
        cwd=args.batch.analysis
    )


def run_bpipe(args):
    check_java(args)
    run(
        f'bpipe {" ".join(args.bpipe_args)}',
        cwd=args.batch.analysis
    )


def main():
    parser = argparse.ArgumentParser(description='Runs core Cpipe functions')
    parser.set_defaults(func=lambda x: parser.print_usage())
    subparsers = parser.add_subparsers(dest='command')

    # Run command
    run_parser = subparsers.add_parser('run', help='Runs the analysis pipeline')
    run_parser.add_argument('batch', type=arg_validation.existing_batch,
                            help='The name of the batch that you want to run, relative to the batches directory')
    run_parser.add_argument('-j', '--no-java-check', action='store_true',
                            help="Don't check that the system Java is a compatible version")
    run_parser.add_argument('bpipe_opts', nargs=argparse.REMAINDER,
                            help="Arguments to pass to the underlying bpipe command")
    run_parser.set_defaults(func=run_pipeline)

    # Batch command
    batch_parser = subparsers.add_parser('batch', help='Manage Cpipe batches and metadata files')
    manage_batch.setup_parser(batch_parser)
    batch_parser.set_defaults(func=manage_batch.execute)

    # Stop command
    stop_parser = subparsers.add_parser('stop', help='Stops a currently running Cpipe batch')
    stop_parser.set_defaults(func=stop_pipeline)

    # Test command
    test_parser = subparsers.add_parser('test', help='Run the cpipe test suite')
    test_parser.add_argument('-j', '--no-java-check', action='store_true',
                             help="Don't check that the system Java is a compatible version")
    test_parser.set_defaults(func=test_pipeline)

    # Design command
    test_parser = subparsers.add_parser('design', help='Create and modify designs and their genelists')
    manage_genelists.setup_parser(test_parser)
    test_parser.set_defaults(func=manage_genelists.execute)

    # Bpipe command
    bpipe_parser = subparsers.add_parser('bpipe', help='Pass a command to bpipe')
    bpipe_parser.add_argument('batch', type=arg_validation.existing_batch,
                              help='The name of the batch, relative to the batches directory, that you want to run bpipe on')
    bpipe_parser.add_argument('-j', '--no-java-check', nargs='?', const=True,
                              help="Don't check that the system Java is a compatible version")
    bpipe_parser.add_argument('bpipe_args', nargs=argparse.REMAINDER,
                              help="Arguments to pass to the underlying bpipe command")
    bpipe_parser.set_defaults(func=run_bpipe)

    # Parse then execute
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
