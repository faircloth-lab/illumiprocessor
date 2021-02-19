#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 31 January 2014 12:26 PST (-0800)
"""


import sys
import configparser
from illumiprocessor import core
from illumiprocessor.log import setup_logging


def main(args):
    core.check_dependencies(args)
    # setup logging
    log, my_name = setup_logging(args)
    # setup config instance
    conf = configparser.ConfigParser()
    # preserve case of entries & read
    conf.optionxform = str
    conf.read(args.config)
    # setup multiprocessing
    pool = core.setup_multiprocessing(args)
    reads = []
    for start_name, end_name in conf.items("names"):
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
        list(map(core.runner, work))
    # push the output down to next line
    print("")
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))
