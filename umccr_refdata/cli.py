import logging
import sys


from . import api
from . import arguments


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)


def entry():
    """Commandline entry point"""
    # Set up logging
    logging.basicConfig()

    # Get command line arguments
    args = arguments.get_arguments()
    LOGGER.debug(f'got arguments {args}')

    # Run required command
    refdata_info = api.read_refdata_information(args.refdata_info_fp)
    if args.subparser_name == 'pull':
        api.pull_bundle(
            args.bundle_name,
            args.output_dir,
            refdata_info,
            args.git_repo_url,
            args.dvc_remote_name,
            args.cache_dir,
            args.git_tag,
        )
    elif args.subparser_name == 'locate':
        filepaths = api.locate_file_paths(
            args.identifier,
            refdata_info,
            args.match_dict,
        )
        if args.reference_data_dir:
            filepaths = api.resolve_refdata_paths(filepaths, args.reference_data_dir)
        print(*filepaths, sep='\n')
    else:
        LOGGER.error(f'Got bad subcommand \'{args.subparser_name}\'')
        sys.exit(1)
