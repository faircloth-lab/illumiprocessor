#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 31 January 2014 11:38 PST (-0800)
"""


import argparse
from illumiprocessor import core
from illumiprocessor.pth import get_user_path

# import pdb


def get_trimmomatic_path():
    return get_user_path("executables", "trimmomatic")


def get_args():
    parser = argparse.ArgumentParser(
        description="Batch trim Illumina reads for adapter contamination and "
        "low quality bases using Trimmomatic",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        required=True,
        help="The input directory of raw reads to trim.",
        action=core.FullPaths,
    )
    parser.add_argument(
        "--output",
        required=True,
        help="The output directory of clean reads to create.",
        action=core.CreateDir,
    )
    parser.add_argument(
        "--config",
        required=True,
        help="A configuration file containing the tag:sample mapping and "
        "renaming options.",
    )
    parser.add_argument(
        "--trimmomatic",
        default=get_trimmomatic_path(),
        action=core.FullPaths,
        help="The path to the trimmomatic-0.XX.jar file.",
    )
    parser.add_argument(
        "--min-len", type=int, default=40, help="The minimum length of reads to keep."
    )
    parser.add_argument(
        "--no-merge",
        action="store_true",
        default=False,
        help="When trimming PE reads, do not merge singleton files.",
    )
    parser.add_argument(
        "--cores", type=int, default=1, help="Number of compute cores to use."
    )
    parser.add_argument(
        "--r1-pattern",
        type=str,
        default=None,
        help="An optional regex pattern to find R1 reads.",
    )
    parser.add_argument(
        "--r2-pattern",
        type=str,
        default=None,
        help="An optional regex pattern to find R2 reads.",
    )
    parser.add_argument(
        "--se",
        action="store_true",
        default=False,
        help="""Single-end reads.""",
    )
    parser.add_argument(
        "--phred",
        type=str,
        choices=("phred33", "phred64"),
        default="phred33",
        help="""The type of fastq encoding used.""",
    )
    parser.add_argument(
        "--log-path",
        action=core.FullPaths,
        type=core.is_dir,
        default=None,
        help="""The path to a directory to hold logs.""",
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use.""",
    )
    return parser.parse_args()


def main():
    from illumiprocessor.main import main

    args = get_args()
    main(args)
