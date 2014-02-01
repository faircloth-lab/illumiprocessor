#!/usr/bin/env python
# encoding: utf-8

from distutils.core import setup

setup(
    name='illumiprocessor',
    version='2.0.6',
    description='Automated Illumina read trimming using trimmomatic',
    url='https://github.com/faircloth-lab/illumiprocessor',
    author='Brant C. Faircloth',
    author_email='borg@faircloth-lab.org',
    license='BSD',
    platforms='any',
    packages=[
        'illumiprocessor',
        'illumiprocessor/cli'
    ],
    scripts=[
        'bin/illumiprocessor'
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    )
