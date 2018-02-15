#!/usr/bin/env python

from setuptools import setup

setup(
    name='haloplexpipe',
    version='0.1',
    author='Jason Steen',
    author_email='jason.steen@monash.edu',
    packages=['src'],
    entry_points={
        'console_scripts': ['haloplexpipe = src.main:main']
    },
    url='https://github.com/SoutheyLab/haloplexpipe',
    license='LICENSE.txt',
    description='hiplexpipe is a bioinformatics pipeline to call variants from HiPlex data.',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.6.3",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
)
