import os
import tempfile
from os.path import join, abspath, dirname, pardir, exists, isdir, isfile
import re
import sys
import yaml
from . import utils


def package_path():
    return dirname(abspath(__file__))


def load_paths():
    """ Depending on the machine name, return a dict conatining system-dependant paths
        to human reference genomes and extras
    """
    with open(join(package_path(), 'paths.yml')) as f:
        loc_by_name = yaml.full_load(f)

    loc_dict = loc_by_name['generic']
    if 'TRAVIS' in os.environ.keys():
        loc_dict = loc_by_name['travis']
    else:
        for ld in loc_by_name.values():
            if 'host_pattern' in ld:
                if re.match(ld['host_pattern'], utils.hostname):
                    loc_dict = utils.update_dict(loc_dict, ld)
                    break
    return loc_dict

loc_dict = load_paths()

name = loc_dict.get('name')
genomes = loc_dict.get('genomes')
genomes_dir = None

# non-reference data stuff
cluster_cmd = loc_dict.get('cluster_cmd')
cluster_jobscript = loc_dict.get('cluster_jobscript')
try:
    ncpus_on_node = loc_dict.get('ncpus_on_node', len(os.sched_getaffinity(0)))
except:
    ncpus_on_node = loc_dict.get('ncpus_on_node', os.cpu_count())


def ref_file_exists(genome='all', key='fa', path=None):
    try:
        return get_ref_file(genome, key, path=path)
    except:
        return False


def get_ref_file(genome='all', key='fa', path=None, must_exist=True):
    """ If path does not exist, checks the "genomes" dictionary for the location.
    """
    if path:
        path = utils.adjust_path(path)
        if exists(path):
            return utils.adjust_path(path)
        else:
            utils.critical(f'Error: {path} is not found as file at {os.getcwd()}')

    g_d = get_genomes_dict(genome)
    d = g_d
    keys = key if not isinstance(key, str) else [key]
    for k in keys:
        d = d.get(k)
        if not d:
            utils.critical(f'Can\'t find refererence file [{".".join(keys)}] for genome "{genome}"'
                     f' for host "{name or utils.hostname}". Available keys: {", ".join(g_d)}')
    if isinstance(d, str):
        path = d
    else:
        utils.critical(f'Error: path in genomes specified by keys [{".".join(keys)}] is not full.'
                 f'Genome: {genome}, cwd: {os.getcwd()}, host: "{name or utils.hostname}"')

    # Resolve found path:
    if not path.startswith('/'):
        if genomes_dir:
            if not isdir(genomes_dir):
                if must_exist:
                    utils.critical(f'Could not find the genomes directory {genomes_dir} for host '
                             f'{name or utils.hostname}')
                return None
            gd = genomes_dir
        else:
            gd = find_genomes_dir()
        path = abspath(join(gd, path))

    path = utils.adjust_path(path)
    if not exists(path):
        if must_exist:
            utils.critical(f'Error: {path} does not exist for genome "{genome}" at host '
                           f'"{name or utils.hostname}"')

    return path


def get_genomes_dict(genome):
    if genome not in genomes:
        utils.critical(f'Error: genome "{genome}" not found for host "{name or utils.hostname}". '
                 f'Available: {", ".join(genomes)}')
    return genomes[genome]


"""
Finding genomes location

A tool may be run with genomes_dir set explicitly via:
    1. --genomes-dir
    2. $UMCCRISE_GENOMES
Or it can have a default for this machine in:
    3. paths.yaml
In these cases the path can be s3 or gds and needs downloading.

If the path is not found anywhere, umccrise may attempt
to find it:
    4. umccrise installation folder under umccrise/genomes
    5. current package installation folder under umccrise/genomes
    6. extras/umccrise/genomes
In these cases the path doesn't need downloading and can be used directly.

Once the path is found and downloaded, .genomes_dir propery is set and can be used directly.
"""


def find_genomes_dir(_new_genomes_url=None):
    """
    Trying:
    1. --genome-dir
    2. $UMCCRISE_GENOMES
    3. paths.yaml
    4. umccrise installation folder under umccrise/genomes
    5. current package installation folder under umccrise/genomes
    6. extras/umccrise/genomes
    """
    global genomes_dir

    tried = []

    if _new_genomes_url:
        # setting up a new genomes location (likely provided by --genomes), even if it was set before
        utils.info(f'Using genomes location that was provided explicitly '
                   f'(in paths.yaml or with --genomes): {_new_genomes_url}')
        genomes_dir = _set_up_new_genomes_dir(_new_genomes_url)
        return genomes_dir
    tried.append(f'--genomes option')

    if genomes_dir is not None:
        # `find_genomes_dir` was called once before, and the genomes dir
        # is already found and/or downloaded locally and can be used directly:
        return genomes_dir

    if os.environ.get('UMCCRISE_GENOMES') is not None:
        gd = os.environ.get('UMCCRISE_GENOMES')
        utils.info(f'Using genomes location set by $UMCCRISE_GENOMES environment variable: {gd}')
        genomes_dir = _set_up_new_genomes_dir(gd)
        return genomes_dir
    tried.append(f'$UMCCRISE_GENOMES env var')

    if 'genomes_dir' in loc_dict:
        gd = loc_dict.get()
        utils.info(f'Using the genomes location from reference_data/paths.yml: {gd}')
        genomes_dir = _set_up_new_genomes_dir(gd)
        return genomes_dir
    tried.append(f'genomes_dir in paths.yaml for location {name or utils.hostname}')

    try:
        from umccrise import package_path as umccrise_path
    except ImportError:
        tried.append(f'umccrsie package parent folder')
    else:
        from umccrise import package_path as umccrise_path
        gd = abspath(join(umccrise_path(), pardir, 'genomes'))
        gd2 = abspath(join(umccrise_path(), pardir, pardir, 'genomes'))
        gd = gd if isdir(gd) else gd2
        if isdir(gd):
            utils.info(f'Using genomes dir at umccrise package parent folder')
            genomes_dir = gd
            return genomes_dir
        tried.append(f'umccrsie package parent folder ({gd})')

    gd = abspath(join(package_path(), pardir, 'genomes'))
    if isdir(gd):
        utils.info(f'Using genomes dir at reference_data package parent folder')
        genomes_dir = gd
        return genomes_dir
    tried.append(f'reference_data package parent folder ({gd})')

    extras = loc_dict.get('extras')
    if extras and isdir(extras):
        gd = join(extras, 'umccrise', 'genomes')
        if isdir(gd):
            utils.info(f'Using genomes dir at production installation of "umccrise" extras')
            genomes_dir = gd
            return genomes_dir
        tried.append(f'production installation of "umccrise" extras ({gd})')

    utils.critical(f'Cannot find "genomes" folder. Tried: {", ".join(tried)}')


def _set_up_new_genomes_dir(genomes_url=None):
    target_genomes_dir = utils.safe_mkdir(utils.adjust_path('~/umccrise_genomes'))
    from . import _version
    version = _version.__version__

    # a local tarball file?
    if isfile(genomes_url) and genomes_url.endswith('.tar.gz') or genomes_url.endswith('.tgz'):
        return utils.extract_tarball_input(genomes_url, target_genomes_dir)

    # a local directory?
    elif isdir(genomes_url):
        return utils.verify_dir(genomes_url, is_critical=True)

    # s3 or gds url?
    else:
        if genomes_url.endswith('/'):  # dropping the trailing /
            genomes_url = genomes_url[:-1]
        versioned_full_genomes_url = f'{genomes_url}_{"".join(version.split("."))}/'
        versioned_genomes_url = f'{genomes_url}_{"".join(version.split(".")[0:2])}/'
        genomes_url = f'{genomes_url}/'  # putting the tralining / back

        urls_to_try = [versioned_full_genomes_url, versioned_genomes_url, genomes_url]

        if genomes_url.startswith('s3://'):
            for url in urls_to_try:
                try:
                    utils.run_simple(f'aws s3 ls {url}')
                except:
                    pass
                else:
                    utils.run_simple(f'aws s3 sync {url} {target_genomes_dir}')
                    return target_genomes_dir
            utils.critical(f'Cannot find reference data at s3. Tried URLs: {urls_to_try}')

        elif genomes_url.startswith('gds://'):
            for url in urls_to_try:
                try:
                    utils.run_simple(f'iap folders list {url} | grep -v -q "No folders found"')
                except:
                    pass
                else:
                    utils.run_simple(f"iap files download {url}/* {target_genomes_dir}")
                    return target_genomes_dir
            utils.critical(f'Cannot find reference data at GDS. Tried URLs: {urls_to_try}')

