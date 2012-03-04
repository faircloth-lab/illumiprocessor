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

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def get_args():
    parser = argparse.ArgumentParser(description='Pre-process Illumina reads')
    parser.add_argument('input', help='The input directory', action=FullPaths)
    parser.add_argument('output', help='The output directory', action=FullPaths)
    parser.add_argument('conf', help='A configuration file containing metadata')
    parser.add_argument('--no-rename',
            dest='rename',
            action='store_false',
            default = True,
            help='Do not rename files using [Map] or [Remap].')
    parser.add_argument('--no-adapter-trim',
            dest='adapter',
            action='store_false',
            default = True,
            help='Do not trim reads for adapter contamination.')
    parser.add_argument('--no-quality-trim', 
            dest='quality',
            action='store_false', 
            default = True, 
            help='Do not trim reads for quality.')
    parser.add_argument('--no-interleave',
            dest='interleave',
            action='store_false',
            default = True, 
            help='Do not interleave trimmed reads.')
    parser.add_argument('--remap', 
            action='store_true',
            default = False,
            help='Remap names onto file using [Remap] section of configuration file.' + \
            ' Used to change file names across many files.')
    parser.add_argument('--se',
            dest='pe',
            action='store_false',
            default = True,
            help='Work with single-end reads (paired-end is default)')
    parser.add_argument('--copy', 
            action='store_true',
            default = False,
            help='Copy, rather than symlink, original files.')
    parser.add_argument('--cleanup',
            action='store_true',
            default = False,
            help='Delete intermediate files.')
    parser.add_argument('--cores', 
            type=int,
            default = 1,
            help='Number of cores to use.')
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

def create_new_dir(base, dirname = None):
    if dirname is None:
        pth = base
    else:
        pth = os.path.join(base, dirname)
    if not os.path.exists(pth):
        os.makedirs(pth)
    return pth


def make_dirs_and_rename_files(inpt, output, sample_map, rename, copy):
    newpths = []
    if rename:
        if not copy:
            print "Symlinking files to output directories...\n"
        else:
            print "Moving files to output directories...\n"
        for old, new in sample_map.iteritems():
            newbase = create_new_dir(output, new)
            newpth = create_new_dir(newbase, 'untrimmed')
            # prep to move files
            oldpth = os.path.join(inpt, '.'.join([old, 'fastq.gz']))
            newpth = os.path.join(newpth, '.'.join([sample_map[old], 'fastq.gz']))
            if not copy:
                os.symlink(oldpth, newpth)
                print "\t{} (sym =>) {}".format(oldpth, newpth)
            else:
                shutil.copyfile(oldpth, newpth)
                print "\t{} => {}".format(oldpth, newpth)
            newpths.append(newbase)
    else:
            for old, new in sample_map.iteritems():
                newbase = os.path.join(output, new)
                newpths.append(newbase)
    return newpths

def build_adapters_file(conf, output):
    adapters = os.path.join(output, 'adapters.fasta')
    if not os.path.exists(adapters):
        f = open(adapters, 'w')
        for adp in conf.items('adapters'):
            f.write(">{}\n{}".format(adp[0], adp[1]))
        f.close()

def scythe_runner(inpt):
    basedir = os.path.split(inpt)[0]
    adapters = os.path.join(basedir, 'adapters.fasta')
    inbase = os.path.basename(inpt)
    infile = ''.join([inbase, ".fastq.gz"])
    inpth = os.path.join(inpt, 'untrimmed', infile)
    outpth = create_new_dir(inpt, 'adapter-trimmed')
    outpth = open(os.path.join(outpth, infile), 'wb')
    statpth = create_new_dir(inpt, 'stats')
    statpth = open(os.path.join(statpth, 'adapter-contam.txt'), 'w')
    
    cmd = ['scythe', '-a', adapters, '-q', 'sanger', inpth]
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
    splitpth = create_new_dir(inpt, 'split-adapter-trimmed')
    infile = ''.join([inbase, ".fastq.gz"])
    inpth = os.path.join(inpt, 'adapter-trimmed', infile)
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
    return

def sickle_pe_runner(inpt):
    inbase = os.path.basename(inpt)
    splitpth = os.path.join(inpt, 'split-adapter-trimmed')
    qualpth = create_new_dir(inpt, 'split-adapter-quality-trimmed')
    # infiles
    r1 = os.path.join(splitpth, ''.join([inbase, "-read1.fastq.gz"]))
    r2 = os.path.join(splitpth, ''.join([inbase, "-read2.fastq.gz"]))
    # outfiles
    out1 = os.path.join(qualpth, ''.join([inbase, "-read1.fastq"]))
    out2 = os.path.join(qualpth, ''.join([inbase, "-read2.fastq"]))
    outS = os.path.join(qualpth, ''.join([inbase, "-read-singleton.fastq"]))
    # make sure we have stat output dir and file
    statpth = create_new_dir(inpt, 'stats')
    statpth = open(os.path.join(statpth, 'sickle-trim.txt'), 'w')
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

def trim_low_qual_reads(pool, newpths, pe = True):
    sys.stdout.write("\nTrimming low quality reads")
    sys.stdout.flush()
    if pe:
        if pool:
            pool.map(sickle_pe_runner, newpths)
        else:
            map(sickle_pe_runner, newpths)
    return

def interleave_reads_runner(inpt):
    inbase = os.path.basename(inpt)
    qualpth = os.path.join(inpt, 'split-adapter-quality-trimmed')
    interpth = create_new_dir(inpt, 'interleaved-adapter-quality-trimmed')

    r1 = os.path.join(qualpth, ''.join([inbase, "-read1.fastq.gz"]))
    r2 = os.path.join(qualpth, ''.join([inbase, "-read2.fastq.gz"]))
    out = os.path.join(interpth, ''.join([inbase, "-read-interleaved.fastq.gz"]))

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
    return

def cleanup_intermediate_files(newpths, interleave):
    dirs = ['adapter-trimmed', 'split-adapter-trimmed']
    for pth in newpths:
        for d in dirs:
            try:
                shutil.rmtree(os.path.join(pth, d))
            except:
                pass
        if interleave:
            shutil.rmtree(os.path.join(pth, 'split-adapter-quality-trimmed'))

def main():
    args = get_args()
    conf = ConfigParser.ConfigParser()
    conf.read(args.conf)
    nproc = multiprocessing.cpu_count()
    if nproc >= 2 and args.cores >= 2:
        pool = multiprocessing.Pool(args.cores)
    else:
        pool = None
    if args.remap:
        names = conf.items('remap')
    else:
        names = conf.items('map')

    create_new_dir(args.output, None)
    sample_map = get_tag_names_from_sample_file(args.input, names)
    newpths = make_dirs_and_rename_files(args.input, args.output, sample_map, args.rename, args.copy)
    if args.adapter:
        build_adapters_file(conf, args.output)
        trim_adapter_sequences(pool, newpths)
    if args.quality:
        if args.pe:
            split_reads(pool, newpths)
            trim_low_qual_reads(pool, newpths)
    if args.pe and args.interleave:
       interleave_reads(pool, newpths)
    if args.cleanup:
        cleanup_intermediate_files(newpths, args.interleave)

    print ""
    
if __name__ == '__main__':
    main()