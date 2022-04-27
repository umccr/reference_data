import argparse
import json
import pathlib
import sys


def get_arguments():
    """Parses and processes command line arguments.

    :returns: Namespace populated with fully processed arguments
    :rtype: argparse.Namespace
    """
    # Parse arguments
    parser = construct_parser()
    args = parser.parse_args()
    # Process arguments
    if not hasattr(args, 'refdata_info_fp') or args.refdata_info_fp is None:
        info_fn = 'data/refdata_information.yaml'
        args.refdata_info_fp = pathlib.Path(__file__).parent / info_fn
    if args.subparser_name == 'locate':
        if args.match_dict is not None:
            args.match_dict = json.loads(args.match_dict)
        if args.reference_data_dir and not args.reference_data_dir.exists():
            parser.error(f'Reference data directory \'{args.reference_data_dir}\' does not exist')
    elif args.subparser_name is None:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return args


def construct_parser():
    """Constructs commandline argument parser.

    :returns: Initialised ArgumentParser instance
    :rtype: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    create_pull_subparser(subparsers)
    create_locate_subparser(subparsers)
    return parser


def create_pull_subparser(subparsers):
    """Create the 'pull' subparser."""
    subparser = subparsers.add_parser(
        'pull',
        description='Pull reference data bundle',
    )

    group_required = subparser.add_argument_group('Required arguments')
    group_required.add_argument(
        '--bundle_name',
        required=True,
        type=str,
        help='Name of reference data bundle to pull',
    )
    group_required.add_argument(
        '--output_dir',
        required=True,
        type=pathlib.Path,
        help='Output directory for DVC repo',
    )

    group_other = subparser.add_argument_group('Other arguments')
    group_other.add_argument(
        '--refdata_info_fp',
        required=False,
        type=pathlib.Path,
        help='Reference data information filepath',
    )
    group_other.add_argument(
        '--git_repo_url',
        required=False,
        type=str,
        help='Git URL for reference data DVC repo',
    )
    group_other.add_argument(
        '--dvc_remote_name',
        required=False,
        type=str,
        help='Name of DVC remote to use',
    )
    group_other.add_argument(
        '--cache_dir',
        required=False,
        type=pathlib.Path,
        help='Local directory to use as DVC cache',
    )
    group_other.add_argument(
        '--git_tag',
        required=False,
        type=str,
        help='Git tag use for checkout prior to DVC pull',
    )


def create_locate_subparser(subparsers):
    """Create the 'locate' subparser."""
    subparser = subparsers.add_parser(
        'locate',
        description='Locate file in local reference data bundle',
    )

    group_required = subparser.add_argument_group('Required arguments')
    group_required.add_argument(
        '--identifier',
        required=True,
        type=str,
        help='File identifier',
    )

    group_other = subparser.add_argument_group('Other arguments')
    group_other.add_argument(
        '--reference_data_dir',
        required=False,
        type=pathlib.Path,
        help='Reference data directory',
    )
    group_other.add_argument(
        '--match_dict',
        required=False,
        type=str,
        help='Key-value pairs to select returned file entries (JSON string)',
    )
    group_other.add_argument(
        '--refdata_info_fp',
        required=False,
        type=pathlib.Path,
        help='Reference data information filepath',
    )
