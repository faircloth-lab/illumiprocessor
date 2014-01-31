#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 31 January 2014 12:43 PST (-0800)
"""

import os
import pytest
import shutil
import ConfigParser
from illumiprocessor import core
from illumiprocessor.cli.main import get_trimmomatic_path


#import pdb


class TruHTFakeArgs:
    def __init__(self):
        self.input = os.path.join(
            os.path.dirname(__file__),
            "truht/raw-reads"
        )
        self.output = os.path.join(
            os.path.dirname(__file__),
            "truht/clean"
        )
        # remove existing clean reads folder
        try:
            shutil.rmtree(self.output)
        except OSError:
            pass
        self.config = os.path.join(
            os.path.dirname(__file__),
            "truht/tru-seq-ht.conf"
        )
        self.trimmomatic = get_trimmomatic_path()
        self.min_len = 40
        self.no_merge = False
        self.cores = 1
        self.se = False
        self.phred = "phred33"
        self.log_path = os.getcwd()
        self.verbosity = "ERROR"
        self.r1_pattern = None
        self.r2_pattern = None


@pytest.fixture(scope="module")
def fake_truht_args(request):
    args = TruHTFakeArgs()
    def clean():
        try:
            shutil.rmtree(args.output)
        except:
            pass
    request.addfinalizer(clean)
    return args


@pytest.fixture(scope="module")
def fake_truht_reads():
    args = TruHTFakeArgs()
    conf = ConfigParser.ConfigParser()
    # preserve case of entries & read
    conf.optionxform = str
    conf.read(args.config)
    reads = []
    for start_name, end_name in conf.items('names'):
        reads.append(core.SequenceData(
            args,
            conf,
            start_name,
            end_name)
        )
    return reads
