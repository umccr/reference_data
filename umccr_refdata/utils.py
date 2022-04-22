import collections
import os
import socket
import subprocess
import time
from os.path import join, abspath, dirname, pardir, exists, isdir, isfile
import sys


def get_hostname():
    return os.environ.get('HOSTNAME') or os.environ.get('HOST') or socket.gethostname()

hostname = get_hostname()


# Expand paths beginning with '~' or '~user'.
# '~' means $HOME; '~user' means that user's home directory.
# If the path doesn't begin with '~', or if the user or $HOME is unknown,
# the path is returned unchanged (leaving error reporting to whatever
# function is called with the expanded path as argument).
# See also module 'glob' for expansion of *, ? and [...] in pathnames.
# (A function should also be defined to do full *sh-style environment
# variable expansion.)
def expanduser(path):
    """Expand ~ and ~user constructs.

    If user or $HOME is unknown, do nothing."""
    if path[:1] != '~':
        return path
    i, n = 1, len(path)
    while i < n and path[i] not in '/\\':
        i = i + 1

    if 'HOME' in os.environ:
        userhome = os.environ['HOME']
    elif 'USERPROFILE' in os.environ:
        userhome = os.environ['USERPROFILE']
    elif not 'HOMEPATH' in os.environ:
        return path
    else:
        try:
            drive = os.environ['HOMEDRIVE']
        except KeyError:
            drive = ''
        userhome = join(drive, os.environ['HOMEPATH'])

    if i != 1:  # ~user
        userhome = join(dirname(userhome), path[1:i])

    return userhome + path[i:]


def remove_quotes(s):
    if s and s[0] in ['"', "'"]:
        s = s[1:]
    if s and s[-1] in ['"', "'"]:
        s = s[:-1]
    return s


def adjust_path(path):
    if path is None: return None

    path = remove_quotes(path)
    if path is None: return None

    path = expanduser(path)
    if path is None: return None

    path = abspath(path)
    if path is None: return None

    return path.replace('//', '/')


def verify_dir(dirpath, description='', silent=False, is_critical=False):
    if dirpath is None:
        msg = (description + ': i' if description else 'I') + 's not specified (None).'
        info(msg, silent, is_critical)
        return dirpath

    if not dirpath:
        msg = (description + ': d' if description else 'D') + 'ir name is empty.'
        info(msg, silent, is_critical)
        return None

    dirpath = adjust_path(dirpath)
    if not exists(dirpath):
        msg = (description + ': ' if description else '') + dirpath + ' does not exist.'
        info(msg, silent, is_critical)
        return None

    if not isdir(dirpath):
        msg = (description + ': ' if description else '') + dirpath + ' is not a directory.'
        info(msg, silent, is_critical)
        return None

    return dirpath


def critical(msg):
    sys.stderr.write(msg + '\n')
    sys.exit(1)


def info(msg, silent=False, is_critical=False):
    if is_critical:
        critical(msg)
    if not silent:
        sys.stderr.write(msg + '\n')


def update_dict(to_update, updates):
    """ Recursively updates a nested dict `to_update` from a nested dict `updates`
    """
    for key, val in updates.items():
        if isinstance(val, collections.Mapping):
            to_update[key] = update_dict(to_update.get(key, {}), val)
        else:
            to_update[key] = updates[key]
    return to_update


def run_simple(cmd, silent=False):
    if not silent:
        info(' '.join(str(x) for x in cmd) if not isinstance(cmd, str) else cmd)
    subprocess.check_call(cmd, shell=True)


def safe_mkdir(dirpath, descriptive_name=''):
    """ Multiprocessing-safely and recursively creates a directory
    """
    if not dirpath:
        critical(f'Path is empty: {descriptive_name if descriptive_name else ""}')

    if isdir(dirpath):
        return dirpath

    if isfile(dirpath):
        critical(descriptive_name + ' ' + dirpath + ' is a file.')

    num_tries = 0
    max_tries = 10

    while not exists(dirpath):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dirpath)
        except OSError as e:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dirpath


def extract_tarball_input(tgz_fpath, target_genomes_dir):
    assert tgz_fpath.endswith('.tar.gz') or tgz_fpath.endswith('.tgz')
    run_simple(f'tar -xzf {tgz_fpath} --directory {target_genomes_dir}')
    assert isdir(target_genomes_dir)
    return target_genomes_dir



