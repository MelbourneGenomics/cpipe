from jinja2 import Template
import os
from cpipe.paths import BASE
from argparse_prompt import PromptParser
import getpass

# Filepaths
in_template = BASE / 'pipeline/bpipe.config.template'
out_template = BASE / 'pipeline/bpipe.config'

def input_with_default(prompt: str, default: str):
    """Asks for user input with the given prompt, and a default value if they don't enter anything"""

    val = input(f'{prompt}\n({default}) ')
    return val or default
            
def main():
    parser = PromptParser()
    parser.add_argument('--executor', '-e', choices=['slurm', 'torque', 'pbspro', 'none'], help='Which queueing system do you use?', default='torque')
    parser.add_argument('--account', '-a', help='What account name do you want to use when running jobs in this system?', default=getpass.getuser())
    parser.add_argument('--queue', '-q', help='Which queue do you want to use when running jobs in this system?', default='main')
    args = parser.parse_args()

    # Write out the new config file
    template = Template(in_template.read_text())
    out_template.write_text(template.render({
        'executor': args.executor,
        'account': args.account,
        'queue': args.queue
    }))

if __name__ == '__main__':
    main()

