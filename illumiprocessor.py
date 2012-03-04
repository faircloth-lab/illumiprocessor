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
 * sickle:     https://github.com/najoshi/sickle
 * seqtools:   https://github.com/faircloth-lab/seqtools

USAGE
------

python illumiprocessor.py \
    indata/ \
    outdata/ \
    pre-process-example.conf

"""

import os
import sys
import shutil
import argparse
import subprocess
import ConfigParser
import multiprocessing

from itertools import izip
from seqtools.sequence import fastq

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Pre-process Illumina reads')
    parser.add_argument('input', help='The input directory')
    parser.add_argument('output', help='The output directory')
    parser.add_argument('conf', help='A configuration file containing metadata')
    parser.add_argument('--remap', 
            action='store_true', default = False,
            help='Remap names onto file using [Remap] section of configuration file.' + \
            ' Used to change file names across many files.')
    parser.add_argument('--move', 
            action='store_true', default = False,
            help='Move, rather than copy, original files.')
    return parser.parse_args()

def get_tag_names_from_sample_file(inpt, names):
    sample_map = {k:v for k,v in names}
    for f in sample_map.keys():
        try:
            fname = '.'.join([f, 'fastq.gz'])
            os.path.isfile(os.path.join(inpt, fname))
        except:
            fname = '.'.join([f, 'fastq'])
            os.path.isfile(os.path.join(inpt, fname))
        finally:
            "File {} does not exist".format(f)
    return sample_map

def make_dirs_and_rename_files(input, output, sample_map, move = False):
    newpths = []
    dirs  = [
            'untrimmed', 
            'adapter-trimmed',
            'split-adapter-trimmed',
            'split-adapter-quality-trimmed',
            'interleaved-adapter-quality-trimmed',
            'stats'
        ]
    if not move:
        print "Copying files to output directories..."
    else:
        print "Moving files to output directories..."
    for old, new in sample_map.iteritems():
        # setup new directories
        newdir = os.path.join(output, new)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        for d in dirs:
            pth = os.path.join(newdir, d)
            if not os.path.exists(pth):
                os.makedirs(pth)
        # prep to move files
        oldpth = os.path.join(input, '.'.join([old, 'fastq.gz']))
        newpth = os.path.join(newdir, 'untrimmed', '.'.join([sample_map[old], 'fastq.gz']))
        if not move:
            shutil.copyfile(oldpth, newpth)
        else:
            shutil.move(oldpth, newpth)
        newpths.append(newdir)
        print "\t{} => {}".format(oldpth, newpth)
    return newpths

def build_adapters_file(conf, output):
    adapters = 'adapters.fasta'
    if not os.path.exists(adapters):
        f = open(adapters, 'w')
        for adp in conf.items('adapters'):
            f.write(">{}\n{}".format(adp[0], adp[1]))
        f.close()

def scythe_runner(inpt):
    inbase = os.path.basename(inpt)
    infile = ''.join([inbase, ".fastq.gz"])
    inpth = os.path.join(inbase, 'untrimmed', infile)
    outpth = open(os.path.join(inbase, 'adapter-trimmed', infile), 'wb')
    statpth = open(os.path.join(inbase, 'stats', 'adapter-contam.txt'), 'w')
    
    cmd = ['scythe', '-a', 'adapters.fasta', '-q', 'sanger', inpth]
    proc1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = statpth)
    proc2 = subprocess.Popen(['gzip'], stdin=proc1.stdout, stdout = outpth)
    proc1.stdout.close()
    output = proc2.communicate()
    outpth.close()
    statpth.close()
    sys.stdout.write(".")
    sys.stdout.flush()

def trim_adapter_sequences(pool, newpths):
    sys.stdout.write("\nTrimming adapter contamination")
    sys.stdout.flush()
    if pool:
        pool.map(scythe_runner, newpths)
    else:
        map(scythe_runner, newpths)
    return

def split_reads_runner(inpt):
    inbase = os.path.basename(inpt)
    splitpth = os.path.join(inbase, 'split-adapter-trimmed')
    if not os.path.exists(splitpth):
        os.makedirs(splitpth)
    infile = ''.join([inbase, ".fastq.gz"])
    inpth = os.path.join(inbase, 'adapter-trimmed', infile)
    reads = fastq.FasterFastqReader(inpth)
    
    out1 = ''.join([inbase, '-read1','.fastq.gz'])
    out2 = ''.join([inbase, '-read2','.fastq.gz'])
    
    r1 = fastq.FasterFastqWriter(os.path.join(splitpth, out1))
    r2 = fastq.FasterFastqWriter(os.path.join(splitpth, out2))

    for read in reads:
        if read[0].split(' ')[1].split(':')[0] == '1':
            r1.write(read)
            first = read[0].split(' ')[0]
        else:
            assert first == read[0].split(' ')[0], "File does not appear interleaved."
            r2.write(read)
    sys.stdout.write('.')
    sys.stdout.flush()
    reads.close()
    r1.close()
    r2.close()

def split_reads(pool, newpths):
    sys.stdout.write("\nSplitting reads")
    sys.stdout.flush()
    if pool:
        pool.map(split_reads_runner, newpths)
    else:
        map(split_reads_runner, newpths)

def sickle_runner(inpt):
    inbase = os.path.basename(inpt)
    splitpth = os.path.join(inbase, 'split-adapter-trimmed')
    qualpth = os.path.join(inbase, 'split-adapter-quality-trimmed')
    # infiles
    r1 = os.path.join(splitpth, ''.join([inbase, "-read1.fastq.gz"]))
    r2 = os.path.join(splitpth, ''.join([inbase, "-read2.fastq.gz"]))
    # outfiles
    out1 = os.path.join(qualpth, ''.join([inbase, "-read1.fastq"]))
    out2 = os.path.join(qualpth, ''.join([inbase, "-read2.fastq"]))
    outS = os.path.join(qualpth, ''.join([inbase, "-read-singleton.fastq"]))
    statpth = open(os.path.join(inbase, 'stats', 'sickle-trim.txt'), 'w')

    #command for sickle (DROPPING ANY Ns)
    cmd = ["sickle", "pe", "-f", r1, "-r", r2,  "-t", "sanger", "-o", out1, "-p", out2, "-s", outS, "-n"]
    proc1 = subprocess.Popen(cmd, stdout=statpth, stderr=subprocess.STDOUT)
    proc1.communicate()
    statpth.close()
    # sickle does not zip, so zip on completion
    for f in [out1, out2, outS]:
        subprocess.Popen(['gzip', f])
    sys.stdout.write(".")
    sys.stdout.flush()

def trim_low_qual_reads(pool, newpths):
    sys.stdout.write("\nTrimming low quality reads")
    sys.stdout.flush()
    if pool:
        pool.map(sickle_runner, newpths)
    else:
        map(sickle_runner, newpths)

def interleave_reads_runner(inpt):
    inbase = os.path.basename(inpt)
    qualpth = os.path.join(inbase, 'split-adapter-quality-trimmed')
    interpth = os.path.join(inbase, 'interleaved-adapter-quality-trimmed')

    r1 = os.path.join(qualpth, ''.join([inbase, "-read1.fastq.gz"]))
    r2 = os.path.join(qualpth, ''.join([inbase, "-read2.fastq.gz"]))
    out = os.path.join(interpth, ''.join([inbase, "read-interleaved.fastq.gz"]))

    read1 = fastq.FasterFastqReader(r1)
    read2 = fastq.FasterFastqReader(r2)
    outfile = fastq.FasterFastqWriter(out)
    for r1,r2 in izip(read1, read2):
        assert r1[0].split(" ")[0] == r2[0].split(" ")[0], \
                "Read FASTQ headers mismatch."
        outfile.write(r1)
        outfile.write(r2)
    sys.stdout.write(".")
    sys.stdout.flush()
    outfile.close()
    read1.close()
    read2.close()

    # move singelton file to interleaved directory
    oldpth = os.path.join(qualpth, ''.join([inbase, "-read-singleton.fastq.gz"]))
    newpth = os.path.join(interpth, ''.join([inbase, "-read-singleton.fastq.gz"]))
    shutil.move(oldpth, newpth)
    # symlink to singleton file from split-adapter-quality-trimmed
    newlink = os.path.join('../', ''.join([inbase, "-read-singleton.fastq.gz"]))
    os.symlink(newlink, oldpth)

def interleave_reads(pool, newpths):
    sys.stdout.write("\nInterleaving reads")
    sys.stdout.flush()
    if pool:
        pool.map(interleave_reads_runner, newpths)
    else:
        map(interleave_reads_runner, newpths)

def main():
    args = get_args()
    conf = ConfigParser.ConfigParser()
    conf.read(args.conf)
    if multiprocessing.nproc > 2:
        pool = multiprocessing.Pool(4)
    else:
        pool = None
    if not args.remap:
        sample_map = get_tag_names_from_sample_file(args.input, conf.items('map'))
    else:
        sample_map = get_tag_names_from_sample_file(args.input, conf.items('remap'))
    newpths = make_dirs_and_rename_files(args.input, args.output, sample_map, args.move)
    # change to working dir
    os.chdir(args.output)
    build_adapters_file(conf, args.output)
    trim_adapter_sequences(pool, newpths)
    split_reads(pool, newpths)
    trim_low_qual_reads(pool, newpths)
    interleave_reads(pool, newpths)
    print ""
    
if __name__ == '__main__':
    main()