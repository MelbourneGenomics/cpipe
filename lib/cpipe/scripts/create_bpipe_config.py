from jinja2 import Template
import os
from cpipe.paths import BASE

# Filepaths
in_template = BASE / 'pipeline/bpipe.config.template'
out_template = BASE / 'pipeline/bpipe.config'

def input_with_default(prompt: str, default: str):
    """Asks for user input with the given prompt, and a default value if they don't enter anything"""

    val = input(f'{prompt}\n({default}) ')
    return val or default
            
def main():
    # Determine queueing system
    executor = input_with_default('What queueing system do you use? "slurm", "torque", "pbspro", or "none" (no queueing system) are accepted.', 'slurm')

    options = {'slurm', "torque", "pbspro", "none"}
    if executor not in options:
        raise Exception(f'Queueing system must be one of: {", ".join(options)}')

    # Get account and queue info
    account = input_with_default('What account name do you want to use when running jobs in this system?', os.getlogin())
    queue = input_with_default('Which queue do you want to use when running jobs in this system?', 'main')

    # Write out the new config file
    template = Template(in_template.read_text())
    out_template.write_text(template.render({
        'executor': executor,
        'account': account,
        'queue': queue
    }))

if __name__ == '__main__':
    main()

