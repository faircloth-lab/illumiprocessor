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
 * scythe:     https://github.com/vsbuffalo/scythe.git
 * sickle:     https://github.com/najoshi/sickle (commit > 09febb6; Feb 28 2012)
 * seqtools:   https://github.com/faircloth-lab/seqtools

USAGE
------

python illumiprocessor.py \
    indata/ \
    outdata/ \
    pre-process-example.conf

"""

import os
import re
import sys
import glob
import errno
import string
import shutil
import logging
import argparse
import subprocess
import ConfigParser
import multiprocessing

import pdb

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def get_args():
    parser = argparse.ArgumentParser(description='Pre-process Illumina reads')
    parser.add_argument(
        '--input',
        required=True,
        help='The input directory',
        action=FullPaths
    )
    parser.add_argument(
        '--output',
        required=True,
        help='The output directory',
        action=FullPaths
    )
    parser.add_argument(
        '--config',
        required=True,
        help='A configuration file containing metadata'
    )
    parser.add_argument(
        '--min-len',
        type=int,
        default=40,
        help='The minimum length of reads to keep (default:40)'
    )
    parser.add_argument(
        '--no-drop-n',
        action='store_false',
        default=True,
        help='Do not drop reads containing ambiguities (default: off)'
    )
    parser.add_argument(
        '--quality-format',
        choices=['sanger', 'solexa', 'illumina'],
        default='sanger',
        help='The quality format of the reads (default: sanger)'
    )
    parser.add_argument(
        '--no-clean',
        action='store_true',
        default=False,
        help='Do not delete intermediate files (default: off)'
    )
    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help='Number of cores to use.'
    )
    return parser.parse_args()


def message():
    print """\n
*****************************************************************
*                                                               *
* Illumiprocessor - automated MPS read trimming                 *
* (c) 2011-2013 Brant C. Faircloth.                             *
* All rights reserved and no guarantees.                        *
*                                                               *
*****************************************************************\n\n"""


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
        self._get_read_data()
        self._get_tag_data(conf)

    def revcomp(self, seq):
        complement = string.maketrans('acgtACGT', 'tgcaTGCA')
        return seq.translate(complement)[::-1]

    def _get_read_data(self):
        all_reads = glob.glob("{}*".format(os.path.join(self.input_dir, self.start_name)))
        for pth in all_reads:
            name = os.path.basename(pth)
            if re.search("{}_(?:.*)_(R1|READ1|Read1|read1)_\d+.fastq(?:.gz)*".format(self.start_name), name):
                self.r1 += (pth,)
            elif re.search("{}_(?:.*)_(R2|READ2|Read2|read2)_\d+.fastq(?:.gz)*".format(self.start_name), name):
                self.r2 += (pth,)

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
            self.i7a = conf.get('adapters', 'i7').replace("*", self.i7s)

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


def build_adapters_file(sample):
    adapters = os.path.join(sample.homedir, 'adapters.fasta')
    with open(adapters, 'w') as outf:
        outf.write(">i5\n{}\n>i7\n{}\n".format(
            sample.i5a,
            sample.i7a
        ))
    return adapters


def scythe_runner(work):
    args, sample = work
    # build sample specific adapters files
    adapters = build_adapters_file(sample)
    # create dir for stat output
    stat_output = os.path.join(sample.homedir, 'stats')
    os.makedirs(stat_output)
    # create dir for program output
    sample.adapt_trimmed = os.path.join(sample.homedir, 'split-adapter-trimmed')
    os.makedirs(sample.adapt_trimmed)
    # get all the read data in the untrimmed dir
    reads = glob.glob(os.path.join(sample.homedir, 'raw-reads', '*.fastq*'))
    for read in reads:
        name = os.path.basename(read)
        with open(os.path.join(sample.adapt_trimmed, name), 'wb') as gzip_file:
            with open(os.path.join(stat_output, '{}-adapter-contam.txt'.format(name)), 'w') as stat_file:
                # Casava >= 1.8 is Sanger encoded - we almost always use Casava >= 1.8
                # set default prior value explicitly, so we know what we used.
                cmd = ["scythe", "-a", adapters, "-q", args.quality_format, "-p", "0.3", read]
                proc1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stat_file)
                proc2 = subprocess.Popen(["gzip"], stdin=proc1.stdout, stdout=gzip_file)
                proc1.stdout.close()
                proc2.communicate()


def sickle_runner(work):
    args, sample = work
    sample.qual_trimmed = os.path.join(sample.homedir, 'split-adapter-quality-trimmed')
    os.makedirs(sample.qual_trimmed)
    # input files
    input = []
    for read in ["READ1", "READ2"]:
        input.append(os.path.join(sample.adapt_trimmed, "{}-{}.fastq.gz".format(
            sample.end_name,
            read
        )))
    # output files
    output = []
    for read in ["READ1", "READ2", "READ-singleton"]:
        output.append(os.path.join(sample.qual_trimmed, "{}-{}.fastq".format(
            sample.end_name,
            read
        )))
    # make sure we have stat output dir and file
    stat_output = os.path.join(sample.homedir, 'stats')
    try:
        os.makedirs(stat_output, True)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(stat_output):
            pass
        else:
            raise
    with open(os.path.join(stat_output, 'sickle-trim.txt'), 'w') as stat_file:
            #command for sickle (DROPPING ANY Ns and seqs < 40 bp)
            cmd = [
                "sickle",
                "pe",
                "-f", input[0],
                "-r", input[1],
                "-t", args.quality_format,
                "-o", output[0],
                "-p", output[1],
                "-s", output[2],
                "-l", str(args.min_len),
            ]
            if not args.no_drop_n:
                cmd += ["-n"]
            proc1 = subprocess.Popen(cmd, stdout=stat_file, stderr=subprocess.STDOUT)
            proc1.communicate()
    # sickle does not gzip, so gzip on completion
    for f in output:
        proc2 = subprocess.Popen(['gzip', f])
        proc2.communicate()


def cleanup_intermediate_files(sample):
    for d in [sample.adapt_trimmed]:
        try:
            shutil.rmtree(d)
        except:
            pass


def runner(work):
    args, sample = work
    # run scythe to trim adapter contamination
    scythe_runner(work)
    sickle_runner(work)
    if not args.no_clean:
        cleanup_intermediate_files(sample)
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
        new_pth = os.path.join(sample.homedir, 'raw-reads')
        os.makedirs(new_pth)
        # link over old data into new dir
        for reads in sample.r1:
            new_file = os.path.join(new_pth, "{}-READ1.fastq.gz".format(sample.end_name))
            os.symlink(reads, new_file)
        for reads in sample.r2:
            new_file = os.path.join(new_pth, "{}-READ2.fastq.gz".format(sample.end_name))
            os.symlink(reads, new_file)

if __name__ == '__main__':
    # display motd
    message()
    # get the arguments from the CLI
    args = get_args()
    # setup config instance
    conf = ConfigParser.ConfigParser()
    # preserve case of entries & read
    conf.optionxform = str
    conf.read(args.config)
    # setup multiprocessing
    pool = setup_multiprocessing(args)
    reads = []
    for start_name, end_name in conf.items('names'):
        reads.append(SequenceData(args, conf, start_name, end_name))
    # create the output directory if not exists
    create_new_dirs(reads)
    # create the set of work for each process
    work = [[args, read] for read in reads]
    # start the cleaning process
    sys.stdout.write("Running")
    sys.stdout.flush()
    if pool:
        pool.map(runner, work)
    else:
        map(runner, work)
