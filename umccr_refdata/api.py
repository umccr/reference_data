import logging
import os
import sys
import yaml


from . import util


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)


# pylint: disable=too-many-arguments,too-many-branches
def pull_bundle(
    bundle_name, output_dir, refdata_info, git_repo_url=None, dvc_remote_name=None, cache_dir=None, git_tag=None
):
    """Pull references data from DVC remote.

    :param str bundle_name: Name of reference data bundle
    :param pathlib.Path output_dir: Output directory
    :param dict refdata_info: Reference data information
    :param str git_repo_url: URL for git repo
    :param str dvc_remote_name: Name of DVC remote
    :param pathlib.Path cache_dir: Cache directory
    :param str git_tag: Git tag for DVC repo
    :returns: None
    :rtype: None
    """
    # Get reference data bundle information
    if 'bundles' not in refdata_info:
        msg = 'Input reference data information is missing the \'bundles\' section'
        LOGGER.error(msg)
        sys.exit(1)
    if not (identifiers := refdata_info['bundles'].get(bundle_name)):
        LOGGER.error(f'Could not find entry for \'{bundle_name}\'')
        sys.exit(1)

    # Get file paths for identifiers
    bundle_filepaths = list()
    for identifier in identifiers:
        filepaths = locate_file_paths(identifier, refdata_info)
        bundle_filepaths.extend(f'reference_data/{fp}' for fp in filepaths)

    # Get some defaults set in reference data information file
    if 'configuration' not in refdata_info:
        msg = 'Input reference data information is missing the \'configuration\' section'
        sys.exit(1)
    configuration = refdata_info['configuration']
    if git_tag is None:
        git_tag = configuration['git_tag']
    if git_repo_url is None:
        git_repo_url = configuration['git_repo_url']
    if dvc_remote_name is None:
        dvc_remote_name = configuration['dvc_remote_name']
    if cache_dir is None:
        cache_dir = output_dir.resolve() / 'cache'

    # Do not checkout with a git tag to point to HEAD
    if git_tag == 'HEAD':
        git_tag = None

    # Clone git repo and checkout specific tag is set
    command_git_clone_base = f'git clone -b {git_tag}' if git_tag else 'git clone'
    command_git_clone = f'{command_git_clone_base} {git_repo_url} {output_dir}'
    util.execute_command(command_git_clone)

    # Set DVC cache directory if specified
    if cache_dir:
        util.execute_command(f'dvc cache dir {cache_dir}', cwd=output_dir)

    # Check remote type
    remote_type = util.get_remote_type(dvc_remote_name, output_dir)
    if remote_type == 'gdrive' and 'GDRIVE_CREDENTIALS_DATA' not in os.environ:
        msg = (
            'Pulling from a gdrive remote currently requires that the \'GDRIVE_CREDENTIALS_DATA\''
            ' environment variable is defined'
        )
        LOGGER.error(msg)
        sys.exit(1)

    # Pull bundle data
    command_dvc_pull = f'dvc pull -r {dvc_remote_name} {" ".join(bundle_filepaths)}'
    util.execute_command(command_dvc_pull, cwd=output_dir)


def locate_file_paths(
    identifier,
    refdata_info,
    match_dict=None,
):
    """Locates reference filepath.

    :param str identifier: Symbolic identifier for a reference file
    :param dict refdata_info: Reference data information
    :param dict match_dict: Optional dict to filter results
    :returns: Filepath
    :rtype: pathlib.Path
    """
    # Search for matching identifier
    if 'files' not in refdata_info:
        msg = f'Input reference data information file \'{refdata_info}\' is missing the \'files\' section'
        LOGGER.error(msg)
        sys.exit(1)
    if not (entries := refdata_info['files'].get(identifier)):
        LOGGER.error(f'Could not find entry for \'{identifier}\'')
        sys.exit(1)

    # Filter results
    if match_dict:
        entries_pass = list()
        for entry in entries:
            # Skip entries where they are missing any match keys or have different values
            match_results = list()
            for k, v in match_dict.items():
                match_results.append(k in entry and entry[k] == v)
            if not all(match_results):
                continue
            entries_pass.append(entry)
        entries = entries_pass
    if len(entries) == 0:
        LOGGER.error(f'Filtered all entries using \'{str(match_dict)}\'')
        sys.exit(1)

    # Enforce selection of a single entry
    if len(entries) > 1:
        msg = f'Found multiple entries for \'{identifier}\''
        if match_dict:
            msg = f'{msg} with match_dict \'{match_dict}\''
        entries_str = '\n\t'.join(str(e) for e in entries)
        msg = f'{msg}:\n\t{entries_str}'
        LOGGER.error(msg)
        sys.exit(1)
    [entry] = entries

    # Set filepaths
    # NOTE(SW): unclear whether there is utility in differeniating 'path' and 'paths'
    path_comp = list()
    if 'path' in entry:
        path_comp.append(entry['path'])
    else:
        path_comp.extend(entry['paths'])
    return path_comp


def read_refdata_information(fp):
    """Load reference file information from disk.

    :param pathlib.Path filepath: Path of reference listing YAML file
    :returns: Dictionary of reference data information
    :rtype: dict
    """
    with fp.open('r') as fh:
        return yaml.safe_load(fh)


def resolve_refdata_paths(paths_comp, refdata_dir):
    """Load reference file information from disk.

    :param list paths_comp: List of reference filepath components
    :param pathlib.Path refdata_dir: Reference data directory
    :returns: Resolved list of reference filepaths
    :rtype: list
    """
    paths = [refdata_dir.resolve() / p for p in paths_comp]
    for path in paths:
        if not path.exists():
            LOGGER.error(f'Filepath \'{path}\' does not exist')
            sys.exit(1)
    return paths
