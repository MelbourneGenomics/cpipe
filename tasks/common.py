from doit.action import CmdAction
import os

HERE = os.path.dirname(__file__)  # The cpipe root directory
ROOT = os.path.dirname(HERE)
TOOLS_ROOT = os.path.join(ROOT, 'tools')
DATA_ROOT = os.path.join(ROOT, 'data')
INSTALL_ROOT = os.path.join(ROOT, 'install')

ENVIRONMENT_FILE = os.path.join(INSTALL_ROOT, 'environment.sh')

def cmd(command, **kwargs):
    """
    Creates a doit CmdAction with certain default parameters used by cpipe, namely using bash as the default shell,
    using the cpipe root as the working directory, and sourcing the environment file before use. Using the environment
    file instead of simply setting the subprocess's env is done to prevent the environment file becoming out of sync
    with the actual environment variables we need set
    :param command: The bash command to run
    :param kwargs: And additional arguments to pass to CmdAction or popen
    :return: A doit CmdAction
    """
    return CmdAction(
        'source {}\n'.format(ENVIRONMENT_FILE) + command,
        executable='bash',
        shell=True,
        cwd=ROOT,
        **kwargs
    )