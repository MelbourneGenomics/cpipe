from manage_batch.parser import create_parser
from manage_batch import list_batches, create_batch, edit_batch, view_batch, validate_metadata, add_sample

parser = create_parser()
args = parser.parse_args()

if args.command == 'list':
    list_batches()
elif args.command == 'create':
    create_batch(args.batch, args.data, args.exome, args.profile, force=args.force, mode=args.mode)
elif args.command == 'edit':
    edit_batch(args.batch, args.editor)
elif args.command == 'view':
    view_batch(args.batch, args.sample)
elif args.command == 'check':
    validate_metadata(args.batch)
elif args.command == 'add_sample':
    add_sample(args.batch, args.data)
else:
    raise ValueError('Unknown command')

