#!/usr/bin/env python
# encoding: utf-8
"""
File: core.py
Author: Brant Faircloth

Created by Brant Faircloth on 30 October 2013 13:10 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

from __future__ import absolute_import
import os
import re
import sys
import glob
import string
import shutil
import hashlib
import argparse
import subprocess
import multiprocessing

import pdb


## TruSeq V2 Kit 2nd set of indexed adapters
iFiveDict= {}
#index 13
iFiveDict["AGTCAA"] = 'CAATCTCGTATGCCGTCTTCTGCTTG'
#index 14
iFiveDict["AGTTCC"] = 'GTATCTCGTATGCCGTCTTCTGCTTG'
#index 15
iFiveDict["ATGTCA"] = 'GAATCTCGTATGCCGTCTTCTGCTTG'
#index 16
iFiveDict["CCGTCC"] = 'CGATCTCGTATGCCGTCTTCTGCTTG'
#index 18
iFiveDict["GTCCGC"] = 'ACATCTCGTATGCCGTCTTCTGCTTG'
#index 19
iFiveDict["GTGAAA"] = 'CGATCTCGTATGCCGTCTTCTGCTTG'
#index 20
iFiveDict["GTGGCC"] = 'TTATCTCGTATGCCGTCTTCTGCTTG'
#index 21
iFiveDict["GTTTCG"] = 'GAATCTCGTATGCCGTCTTCTGCTTG'
#index 22
iFiveDict["CGTACG"] = 'TAATCTCGTATGCCGTCTTCTGCTTG'
#index 23
iFiveDict["GAGTGG"] = 'ATATCTCGTATGCCGTCTTCTGCTTG'
#index 25
iFiveDict["ACTGAT"] = 'ATATCTCGTATGCCGTCTTCTGCTTG'
#index 27
iFiveDict["ATTCCT"] = 'TTATCTCGTATGCCGTCTTCTGCTTG'




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


class SequenceData():
    def __init__(self, args, conf, start_name, end_name):
        self.input_dir = args.input
        self.output_dir = args.output
        self.start_name = start_name
        self.end_name = end_name
        self.homedir = os.path.join(args.output, end_name)
        self.r1 = ()
        self.r2 = ()
        self.i5 = None
        self.i7 = None
        self.i5s = None
        self.i5s_revcomp = False
        self.i7s = None
        self.i5a = None
        self.i7a = None
        self.se = args.se
        if args.r1_pattern is None:
            self.r1_pattern = "{}_(?:.*)_(R1|READ1|Read1|read1)_\d+.fastq(?:.gz)*"
        else:
            self.r1_pattern = args.r1_pattern
        if args.r2_pattern is None:
            self.r2_pattern = "{}_(?:.*)_(R2|READ2|Read2|read2)_\d+.fastq(?:.gz)*"
        else:
            self.r2_pattern = args.r2_pattern
        self._get_read_data()
        self._get_tag_data(conf)

    def revcomp(self, seq):
        complement = string.maketrans('acgtACGT', 'tgcaTGCA')
        return seq.translate(complement)[::-1]

    def _get_read_data(self):
        all_reads = glob.glob("{}*".format(os.path.join(
            self.input_dir, self.start_name)))
        for pth in all_reads:
            name = os.path.basename(pth)
            if re.search(self.r1_pattern.format(self.start_name), name):
                self.r1 += (pth,)
            elif re.search(self.r2_pattern.format(self.start_name), name):
                self.r2 += (pth,)
        if not self.se and (self.r1 == () or self.r2 == ()):
            raise IOError("There is a problem with the read names for {}. "
                          "Ensure you do not have spelling/capitalization "
                          "errors in your conf file.".format(self.start_name))
        elif self.se and self.r1 == ():
            raise IOError("There is a problem with the read names for {}. "
                          "Ensure you do not have spelling/capitalization "
                          "errors in your conf file.".format(self.start_name))

    def _get_tag_data(self, conf):
        tags = dict(conf.items('tag sequences'))
        tag_map = dict(conf.items('tag map'))
        # get the tags for the sample
        combo = tag_map[self.start_name]
        if "," in combo:
            tag1, tag2 = combo.split(',')
            for t in [tag1, tag2]:
                if 'i5' in t:
                    self.i5 = t
                elif 'i7' in t:
                    self.i7 = t
            self.i5s = self.revcomp(tags[self.i5])
            self.i5s_revcomp = True
            self.i5a = conf.get('adapters', 'i5').replace("*", self.i5s)
            self.i7s = tags[self.i7]
            self.i7a = conf.get('adapters', 'i7').replace("*", self.i7s)
        else:
            self.i7 = combo
            self.i7s = tags[self.i7]
            p1, p2 = conf.get('adapters', 'i7').split('*')
            if self.i7s in iFiveDict:
                self.i7a = p1 + self.i7s + iFiveDict[self.i7s]
            else:
                self.i7a = conf.get('adapters', 'i7').replace("*", self.i7s)
            # there is no index sequence in this adapter
            if not self.se:
                self.i5a = conf.get('adapters', 'i5')

    def __repr__(self):
        return "<{}.{} object at {}, name={}; i7={},{}; i5={},{}; i5revcomp={}>".format(
            self.__module__,
            self.__class__.__name__,
            hex(id(self)),
            self.start_name,
            self.i7,
            self.i7s,
            self.i5,
            self.i5s,
            self.i5s_revcomp
        )


def build_se_adapters_file(sample):
    adapters = os.path.join(sample.homedir, 'adapters.fasta')
    with open(adapters, 'w') as outf:
        outf.write(">i7\n{}\n".format(
            sample.i7a
        ))
    return adapters


def build_pe_adapters_file(sample):
    adapters = os.path.join(sample.homedir, 'adapters.fasta')
    with open(adapters, 'w') as outf:
        outf.write(">i5\n{}\n>i7\n{}\n".format(
            sample.i5a,
            sample.i7a
        ))
    return adapters


def trimmomatic_se_runner(work):
    args, sample = work
    # build sample specific adapters files
    adapters = build_se_adapters_file(sample)
    # create dir for stat output
    stat_output = os.path.join(sample.homedir, 'stats')
    os.makedirs(stat_output)
    # create dir for program output
    sample.trimmed = os.path.join(
        sample.homedir,
        'split-adapter-quality-trimmed'
    )
    os.makedirs(sample.trimmed)
    # get all the read data in the untrimmed dir
    input = []
    for read in ["READ1"]:
        input.append(os.path.join(sample.raw_reads, "{}-{}.fastq.gz".format(
            sample.end_name,
            read
        )))
    output = []
    for read in ["READ1"]:
        output.append(os.path.join(sample.trimmed, "{}-{}.fastq.gz".format(
            sample.end_name,
            read
        )))
    with open(os.path.join(stat_output, '{}-adapter-contam.txt'.format(sample.end_name)), 'w') as stat_file:
        # Casava >= 1.8 is Sanger encoded "-phred33" - we almost always use
        # Casava >= 1.8
        cmd = [
            "java",
            "-jar",
            args.trimmomatic,
            "SE",
            "-{}".format(args.phred),
            input[0],
            output[0],
            "ILLUMINACLIP:{}:2:30:10".format(adapters),
            "LEADING:5",
            "TRAILING:15",
            "SLIDINGWINDOW:4:15",
            "MINLEN:{}".format(args.min_len)
        ]
        proc1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stat_file)
        proc1.communicate()


class TestJavaAndTrimmomatic:
    def __init__(self):
        cmd = ["java", "-version"]
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        self.stdout, self.stderr = proc.communicate()
        self.result_split = self.stderr.strip().split("\n")

    def get_version(self, s):
        result = re.search('java version "(.*)"', s)
        return result.groups()[0].split(".")

    def test_install(self):
        assert self.result_split[0].startswith("java version"), \
            "Java does not appear to be installed"

    def test_version(self):
        version = self.get_version(self.result_split[0])
        assert (version[:2] == ["1", "6"] or version[:2] == ["1", "7"]), \
            "Java does not appear to be 1.6 or 1.7"

    def test_trimmomatic(self, trimmomatic_pth):
        cmd = ["java", "-jar", trimmomatic_pth]
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        self.stdout, self.stderr = proc.communicate()
        # Looks like trimmomatic does not report a version number in any
        # read way, so just check start of output to ensure it appears to be
        # what we want.
        assert self.stderr.startswith("Usage: "), \
            "Trimmomatic does not appear to be installed"


def check_dependencies(args):
    java = TestJavaAndTrimmomatic()
    java.test_install()
    java.test_version()
    java.test_trimmomatic(args.trimmomatic)


def trimmomatic_pe_runner(work):
    args, sample = work
    # build sample specific adapters files
    adapters = build_pe_adapters_file(sample)
    # create dir for stat output
    stat_output = os.path.join(sample.homedir, 'stats')
    os.makedirs(stat_output)
    # create dir for program output
    sample.trimmed = os.path.join(
        sample.homedir,
        'split-adapter-quality-trimmed'
    )
    os.makedirs(sample.trimmed)
    # get all the read data in the untrimmed dir
    input = []
    for read in ["READ1", "READ2"]:
        input.append(os.path.join(sample.raw_reads, "{}-{}.fastq.gz".format(
            sample.end_name,
            read
        )))
    output = []
    for read in ["READ1", "READ1-single", "READ2", "READ2-single"]:
        output.append(os.path.join(sample.trimmed, "{}-{}.fastq.gz".format(
            sample.end_name,
            read
        )))
    with open(os.path.join(stat_output, '{}-adapter-contam.txt'.format(sample.end_name)), 'w') as stat_file:
        # Casava >= 1.8 is Sanger encoded "-phred33" - we almost always use
        # Casava >= 1.8
        cmd = [
            "java",
            "-jar",
            args.trimmomatic,
            "PE",
            "-{}".format(args.phred),
            input[0],
            input[1],
            output[0],
            output[1],
            output[2],
            output[3],
            "ILLUMINACLIP:{}:2:30:10".format(adapters),
            "LEADING:5",
            "TRAILING:15",
            "SLIDINGWINDOW:4:15",
            "MINLEN:{}".format(args.min_len)
        ]
        proc1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stat_file)
        proc1.communicate()


def trimmomatic_merger(sample):
    # rename READ1-single to READ-singleton
    singles = []
    for read in ["READ1-single", "READ2-single"]:
        singles.append(os.path.join(sample.trimmed, "{}-{}.fastq.gz".format(
            sample.end_name,
            read
        )))
    # cat contents of READ2-single into READ-singleton
    new_name = "{}-READ-singleton.fastq.gz".format(sample.end_name)
    new_pth = os.path.join(sample.trimmed, new_name)
    os.rename(singles[0], new_pth)
    # remove READ2-single
    with open(new_pth, 'ab') as outfile:
        shutil.copyfileobj(open(singles[1]), outfile)
    os.remove(singles[1])


def runner(work):
    args, sample = work
    if not args.se:
        # run trimmomatic to trim adapter contamination and low qual bases
        trimmomatic_pe_runner(work)
        if not args.no_merge:
            trimmomatic_merger(sample)
    else:
        trimmomatic_se_runner(work)
    sys.stdout.write(".")
    sys.stdout.flush()


def setup_multiprocessing(args):
    nproc = multiprocessing.cpu_count()
    if nproc >= 2 and args.cores >= 2:
        pool = multiprocessing.Pool(args.cores)
    else:
        pool = None
    return pool


def create_new_dirs(reads):
    # first create the basic output directory - just
    # use first sequence object to get output dir
    if not os.path.exists(reads[0].output_dir):
        os.makedirs(reads[0].output_dir)
    # now create the directory to hold our links to old data
    for sample in reads:
        sample.raw_reads = os.path.join(sample.homedir, 'raw-reads')
        os.makedirs(sample.raw_reads)
        # link over old data into new dir
        for reads in sample.r1:
            new_file = os.path.join(sample.raw_reads, "{}-READ1.fastq.gz".format(sample.end_name))
            os.symlink(reads, new_file)
        for reads in sample.r2:
            new_file = os.path.join(sample.raw_reads, "{}-READ2.fastq.gz".format(sample.end_name))
            os.symlink(reads, new_file)
