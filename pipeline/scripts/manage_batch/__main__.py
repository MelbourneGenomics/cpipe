#!/usr/bin/env python3
import sys

from manage_batch.parser import create_parser
from cpipe_util import Batch

parser = create_parser()
args = parser.parse_args()

if args.command == 'list':
    for batch in Batch.list_all():
        print(batch.name)
elif args.command == 'create':
    Batch.create(args.batch, args.data, args.exome, args.profile, force=args.force, mode=args.mode)
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

