from manage_batch.parser import create_parser

parser = create_parser()
args = parser.parse_args()
args.func(args)
