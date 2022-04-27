import logging
import pathlib
import subprocess
import sys
import urllib.parse


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)


def execute_command(command, cwd=None):
    """Executes commands using subprocess and checks return code.

    :param str command: Command to execute
    :param pathlib.Path cwd: Path to change into prior to execution
    :returns: Data of executed command, including standard streams
    :rtype: subprocess.CompletedProcess
    """
    # pylint: disable=subprocess-run-check
    LOGGER.debug(f'command to execute: {command}')
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        encoding='utf-8',
        cwd=cwd,
    )
    if result.returncode != 0:
        message = (
            f'Got non-zero exit code for: {result.args}\n' f'stdout: {result.stdout}\n' f'stderr: {result.stderr}\n'
        )
        LOGGER.critical(message)
        sys.exit(1)
    LOGGER.debug(f'successfully executed command: {command}')
    return result


def get_remote_type(name, dvc_dp):
    """For a given DVC remote, obtain remote type.

    :param str name: DVC remote name
    :param pathlib.Path dvc_dp: Directory path of the DVC repo
    :returns: DVC remote type
    :rtype: str
    """
    # pylint: disable=useless-else-on-loop
    result = execute_command('dvc remote list', cwd=dvc_dp)
    for line in result.stdout.rstrip().splitlines():
        remote_name, url = line.split('\t')
        if remote_name != name:
            continue
        parsed_url = urllib.parse.urlparse(url)
        return parsed_url.scheme if parsed_url.scheme != '' else 'local'
    else:
        LOGGER.error(f'Did not find DVC remote name \'{name}\' in output:\n{result.stdout}')
        sys.exit(1)


def get_refdata_information_fp():
    """Get the filepath of the reference data information YAML file.

    :returns: Reference data information filepath
    :rtype: pathlib.Path
    """
    info_fn = 'data/refdata_information.yaml'
    return pathlib.Path(__file__).parent / info_fn
