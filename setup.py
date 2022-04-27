#!/usr/bin/env python
# pylint: disable=consider-using-with,unspecified-encoding


import setuptools


import umccr_refdata


setuptools.setup(
    name='umccr_refdata',
    version=umccr_refdata.__version__,
    description='UMCCR reference data API',
    long_description=open('README.md', 'r').read(),
    long_description_content_type='text/markdown',
    author='UMCCR and Contributors',
    author_email='services@umccr.org',
    entry_points={
        'console_scripts': ['umccr_refdata=umccr_refdata.cli:entry'],
    },
    packages=setuptools.find_packages(),
    include_package_data=True,
    python_requires='>=3.8',
    license='GPL',
    url='https://github.com/umccr/reference_data',
)
