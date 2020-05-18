#!/usr/bin/env python
from setuptools import setup

import reference_data
pkg = reference_data.__name__

import versionpy
version = versionpy.get_version(pkg)

setup(
    name=pkg,
    version=version,
    author='Vlad Saveliev',
    author_email='vladislav.saveliev@unimelb.edu.au',
    description='Versioning of reference data used in UMCCR pipelines, and python API to access it',
    keywords='bioinformatics',
    url='https://github.com/umccr/' + pkg,
    license='GPLv3',
    packages=[
        pkg,
    ],
    package_data={
        pkg: versionpy.find_package_files('', pkg),
    },
    include_package_data=True,
    zip_safe=False,
    install_requires=versionpy.get_reqs(),
)
