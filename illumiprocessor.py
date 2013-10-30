#!/usr/bin/env python
# encoding: utf-8

"""
File: illumiprocessor.py
Author: Brant Faircloth

Created by Brant Faircloth on 26 May 2011 14:03 PST (-0800)
Copyright (c) 2011-2012 Brant C. Faircloth. All rights reserved.

Description:


REQUIRES
--------

 * python 2.7
 * Trimmomatic

"""

import os
import sys
import shutil
import argparse
import ConfigParser
from illumiprocessor import core
from illumiprocessor.log import setup_logging

#import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class CreateDir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # get the full path
        d = os.path.abspath(os.path.expanduser(values))
        # check to see if directory exists
        if os.path.exists(d):
            answer = raw_input("[WARNING] Output directory exists, REMOVE [Y/n]? ")
            if answer == "Y":
                shutil.rmtree(d)
            else:
                print "[QUIT]"
                sys.exit()
        # create the new directory
        os.makedirs(d)
        # return the full path
        setattr(namespace, self.dest, d)


def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def get_args():
    parser = argparse.ArgumentParser(
        description='Batch trim Illumina reads for adapter contamination and low quality bases using Trimmomatic',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--input',
        required=True,
        help='The input directory of raw reads to trim.',
        action=FullPaths
    )
    parser.add_argument(
        '--output',
        required=True,
        help='The output directory of clean reads to create.',
        action=CreateDir
    )
    parser.add_argument(
        '--config',
        required=True,
        help='A configuration file containing the tag:sample mapping and renaming options.'
    )
    parser.add_argument(
        '--trimmomatic',
        required=True,
        default="~/bin/trimmomatic-0.30.jar",
        action=FullPaths,
        help='The path to the trimmomatic-0.XX.jar file.'
    )
    parser.add_argument(
        '--min-len',
        type=int,
        default=40,
        help='The minimum length of reads to keep.'
    )
    parser.add_argument(
        '--no-merge',
        action='store_true',
        default=False,
        help='When trimming PE reads, do not merge singleton files.'
    )
    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help='Number of compute cores to use.'
    )
    parser.add_argument(
        '--r1-pattern',
        type=str,
        default=None,
        help='An optional regex pattern to find R1 reads.'
    )
    parser.add_argument(
        '--r2-pattern',
        type=str,
        default=None,
        help='An optional regex pattern to find R2 reads.'
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
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use."""
    )
    return parser.parse_args()


if __name__ == '__main__':
    # get the arguments from the CLI
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    # setup config instance
    conf = ConfigParser.ConfigParser()
    # preserve case of entries & read
    conf.optionxform = str
    conf.read(args.config)
    # setup multiprocessing
    pool = core.setup_multiprocessing(args)
    reads = []
    for start_name, end_name in conf.items('names'):
        reads.append(core.SequenceData(args, conf, start_name, end_name))
    # create the output directory if not exists
    core.create_new_dirs(reads)
    # create the set of work for each process
    work = [[args, read] for read in reads]
    # let us know that read cleaning has started
    log.info("Trimming samples with Trimmomatic")
    # start the cleaning process and output dots to indicate file progress
    sys.stdout.write("Running")
    sys.stdout.flush()
    if pool is not None:
        pool.map(core.runner, work)
    else:
        map(core.runner, work)
    # push the output down to next line
    print ""
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))
