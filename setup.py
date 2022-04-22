#!/usr/bin/env python
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
    packages=setuptools.find_packages(),
    package_data={'umccr_refdata': ['paths.yml']},
    python_requires='>=3.8',
    license='GPL',
    url='https://github.com/umccr/reference_data',
)
