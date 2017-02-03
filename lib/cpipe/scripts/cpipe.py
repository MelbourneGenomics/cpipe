#!/usr/bin/env python3
import argparse
from subprocess import check_output
from os import environ
from cpipe import arg_validation
from collections.abc import Sequence
from doit.doit_cmd import DoitMain
from cpipe.scripts import manage_batch
from cpipe.scripts import manage_genelists
from cpipe import paths


def run(cmd, **kwargs):
    if isinstance(cmd, Sequence) and not isinstance(cmd, str):
        check_output(cmd, cwd=environ.get('CPIPE_ROOT'), **kwargs)
    else:
        check_output(cmd, shell=True, cwd=environ.get('CPIPE_ROOT'), **kwargs)


def check_java(args):
    """Run the java check if necessary"""
    if args.check_java:
        DoitMain().run(['check_java'])


def run_pipeline(args):
    check_java(args)
    args.batch.analysis.mkdir(exist_ok=True, parents=True)
    run(
        f'${paths.BASE}/bpipe run ${args.bpipe_opts} ${paths.BASE}/pipeline/pipeline.groovy {args.batch.metadata_file}',
        cwd=args.batch.analysis
    )


def test_pipeline(args):
    check_java(args)
    run(paths / 'pipeline/scripts/run_unit_tests.sh')
    run(f"{paths / 'pipeline/scripts/run_tests.sh'} detect_mutations_test")


def main():
    parser = argparse.ArgumentParser(description='Runs core Cpipe functions')
    subparsers = parser.add_subparsers(dest='command')

    # Run command
    run_parser = subparsers.add_parser('run', help='Runs the analysis pipeline')
    run_parser.add_argument('batch', type=arg_validation.existing_batch)
    run_parser.add_argument('--no-java-check', '-j', type=arg_validation.existing_batch)
    run_parser.add_argument('bpipe_opts', type=arg_validation.existing_batch, nargs=argparse.REMAINDER)
    run_parser.set_defaults(func=run_pipeline)

    # Batch command
    batch_parser = subparsers.add_parser('batch', help='Manage Cpipe batches and metadata files')
    manage_batch.setup_parser(batch_parser)
    batch_parser.set_defaults(func=manage_batch.execute)

    # Test command
    test_parser = subparsers.add_parser('test', help='Run the cpipe test suite')
    test_parser.set_defaults(func=test_pipeline)

    # Design command
    test_parser = subparsers.add_parser('designs', help='Create and modify designs and their genelists')
    manage_genelists.setup_parser(test_parser)
    test_parser.set_defaults(func=manage_genelists.execute)

    # Parse then execute
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
