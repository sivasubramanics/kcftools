#!/usr/bin/env python
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

with open('README.md') as f:
    long_description = f.read()

setup(
    name='kcftools',
    version='0.1',
    description='python package for IBS and variable k-mer analysis',
    author='Sivasubramani Selvanayagam',
    author_email='siva.selvanayagam@wur.nl',
    # packages=['kcftools'],
    packages=find_packages(),
    install_requires=required,
    long_description=long_description,
    long_description_content_type='text/markdown',
    entry_points={
        'console_scripts': [
            'kcftools = kcftools.main:main'
        ]
    },
)