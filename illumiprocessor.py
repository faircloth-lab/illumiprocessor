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
            default=True,
            help='Do not rename files using [Map] or [Remap].'
        )
    parser.add_argument('--complex',
            dest='complex',
            action='store_true',
            default=False,
            help='Complex file names or two-reads for PE data.'
        )
    parser.add_argument('--no-adapter-trim',
            dest='adapter',
            action='store_false',
            default=True,
            help='Do not trim reads for adapter contamination.'
        )
    parser.add_argument('--no-quality-trim',
            dest='quality',
            action='store_false',
            default=True,
            help='Do not trim reads for quality.'
        )
    parser.add_argument('--no-interleave',
            dest='interleave',
            action='store_false',
            default=True,
            help='Do not interleave trimmed reads.'
        )
    parser.add_argument('--remap',
            action='store_true',
            default=False,
            help='Remap names onto file using [remap] section of configuration file.' + \
            ' Used to change file names across many files.'
        )
    parser.add_argument('--se',
            dest='pe',
            action='store_false',
            default=True,
            help='Work with single-end reads (paired-end is default)'
        )
    parser.add_argument('--copy',
            action='store_true',
            default=False,
            help='Copy, rather than symlink, original files.'
        )
    parser.add_argument('--cleanup',
            action='store_true',
            default=False,
            help='Delete intermediate files.'
        )
    parser.add_argument('--only-cleanup',
            action='store_true',
            default=False,
            help='Delete intermediate files.'
        )
    parser.add_argument('--cores',
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
* (c) 2011-2012 Brant C. Faircloth.                             *
* All rights reserved and no guarantees.                        *
*                                                               *
*****************************************************************\n\n"""


def build_file_name(f, opts, read, directory):
    # use glob here for wilcard expansion
    pth = glob.glob(os.path.join(directory, read.format(name=f)))
    # make sure we get back only 1 match
    assert len(pth) == 1, "Your name format matches more than one file"
    return pth[0]


def check_read_names(opts, f, inpt):
    if opts.tworeads:
        r = [opts.read1, opts.read2]
    else:
        r = [opts.read1]
    for i in r:
        fname = build_file_name(f, opts, i, inpt)
        if not os.path.isfile(fname):
            msg = "{} does not exist".format(fname)
            raise IOError(msg)


def get_tag_names_from_sample_file(inpt, names, remap, opts):
    if remap:
        sample_map = {k:v for k, v in names}
    else:
        sample_map = {k:k for k, v in names}
    for f in sample_map.keys():
        check_read_names(opts, f, inpt)
    return sample_map


def create_new_dir(base, dirname=None):
    if dirname is None:
        pth = base
    else:
        pth = os.path.join(base, dirname)
    if not os.path.exists(pth):
        os.makedirs(pth)
    return pth


def make_dirs_and_rename_files(inpt, output, sample_map, rename, copy, opts):
    newpths = []
    if opts.tworeads:
        reads = [opts.read1, opts.read2]
        regex = re.compile('._(?:R|Read|READ)(\d)_{0,1}.')
    else:
        reads = [opts.read1]
    if rename:
        if not copy:
            print "Symlinking files to output directories...\n"
        else:
            print "Moving files to output directories...\n"
        for old, new in sample_map.iteritems():
            newbase = create_new_dir(output, new)
            newpth = create_new_dir(newbase, 'untrimmed')
            for read in reads:
                # prep to move files
                oldfile = build_file_name(old, opts, read, inpt)
                #pdb.set_trace()
                if len(reads) == 1:
                    newfile = os.path.join(newpth, '.'.join([sample_map[old], 'fastq.gz']))
                else:
                    rnum = regex.search(read).groups()[0]
                    name = "{}-READ{}".format(sample_map[old], rnum)
                    newfile = os.path.join(newpth, '.'.join([name, 'fastq.gz']))
                if not copy:
                    os.symlink(oldfile, newfile)
                    print "\t{} (sym =>) {}".format(oldfile, newfile)
                else:
                    shutil.copyfile(oldfile, newfile)
                    print "\t{} => {}".format(oldfile, newfile)
            newpths.append(newbase)
    else:
        for old, new in sample_map.iteritems():
            newbase = os.path.join(output, new)
            newpths.append(newbase)
    return newpths


def build_adapters_file(conf, inpt):
    #pdb.set_trace()
    outfile = os.path.join(inpt, 'adapters.fasta')
    name = os.path.basename(inpt)
    combos = [i for i in conf.get('combos', name).split(',')]
    seqs = dict(conf.items('indexes'))
    adapters = dict(conf.items('adapters'))
    indexed = {}
    for combo in combos:
        if combo.startswith('n7'):
            indexed[combo] = adapters['n7'].replace('*', seqs[combo])
        elif combo.startswith('n5'):
            indexed[combo] = adapters['n5'].replace('*', seqs[combo])
        elif 'truseq' in combo or 'idt' in combo:
            for k, v in adapters.iteritems():
                indexed[k] = v.replace('*', seqs[combo])
            break
        else:
            for k, v in adapters.iteritems():
                indexed[k] = v.replace('*', seqs[combo])
            break
    if not indexed:
        indexed = adapters
    if not os.path.exists(outfile):
        f = open(outfile, 'w')
        for k, v in indexed.iteritems():
            f.write(">{}\n{}\n".format(k, v))
        f.close()
    return outfile


def scythe_runner(inpt):
    #pdb.set_trace()
    inpt, conf, outdirname = inpt
    # build sample specific adapters files
    adapters = build_adapters_file(conf, inpt)
    topdir = os.path.split(inpt)[0]
    inbase = os.path.basename(inpt)
    infiles = glob.glob(os.path.join(inpt, 'untrimmed', '*.fastq*'))
    for infile in infiles:
        infile_name = os.path.basename(infile)
        outpth = create_new_dir(inpt, outdirname)
        outpth = open(os.path.join(outpth, infile_name), 'wb')
        statpth = create_new_dir(inpt, 'stats')
        statpth = open(
                os.path.join(statpth, '{}-adapter-contam.txt'.format(infile_name)), 'w'
            )
        cmd = ['scythe', '-a', adapters, '-q', 'sanger', infile]
        proc1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=statpth)
        proc2 = subprocess.Popen(['gzip'], stdin=proc1.stdout, stdout=outpth)
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
    out1 = ''.join([inbase, '-READ1', '.fastq.gz'])
    out2 = ''.join([inbase, '-READ2', '.fastq.gz'])
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
    #pdb.set_trace()
    inname = os.path.split(inpt)[1]
    splitpth = os.path.join(inpt, 'split-adapter-trimmed')
    qualpth = create_new_dir(inpt, 'split-adapter-quality-trimmed')
    # infiles
    r1 = os.path.join(splitpth, ''.join([inname, "-READ1.fastq.gz"]))
    r2 = os.path.join(splitpth, ''.join([inname, "-READ2.fastq.gz"]))
    # outfiles
    out1 = os.path.join(qualpth, ''.join([inname, "-READ1.fastq"]))
    out2 = os.path.join(qualpth, ''.join([inname, "-READ2.fastq"]))
    outS = os.path.join(qualpth, ''.join([inname, "-READ-singleton.fastq"]))
    # make sure we have stat output dir and file
    statpth = create_new_dir(inpt, 'stats')
    statpth = open(os.path.join(statpth, 'sickle-trim.txt'), 'w')
    #command for sickle (DROPPING ANY Ns)
    cmd = ["sickle", "pe", "-f", r1, "-r", r2,  "-t", "sanger", "-o", out1, "-p", out2, "-s", outS, "-n"]
    proc1 = subprocess.Popen(cmd, stdout=statpth, stderr=subprocess.STDOUT)
    err = proc1.communicate()
    statpth.close()
    # sickle does not zip, so zip on completion
    for f in [out1, out2, outS]:
        proc2 = subprocess.Popen(['gzip', f])
        proc2.communicate()
    sys.stdout.write(".")
    sys.stdout.flush()


def trim_low_qual_reads(pool, newpths, pe=True):
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

    r1 = os.path.join(qualpth, ''.join([inbase, "-READ1.fastq.gz"]))
    r2 = os.path.join(qualpth, ''.join([inbase, "-READ2.fastq.gz"]))
    out = os.path.join(interpth, ''.join([inbase, "-READ1and2-interleaved.fastq.gz"]))

    read1 = fastq.FasterFastqReader(r1)
    read2 = fastq.FasterFastqReader(r2)
    outfile = fastq.FasterFastqWriter(out)
    for r1, r2 in izip(read1, read2):
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
    oldpth = os.path.join(qualpth, ''.join([inbase, "-READ-singleton.fastq.gz"]))
    newpth = os.path.join(interpth, ''.join([inbase, "-READ-singleton.fastq.gz"]))
    shutil.move(oldpth, newpth)
    # symlink to singleton file from split-adapter-quality-trimmed
    newlink = os.path.join('../interleaved-adapter-quality-trimmed', ''.join([inbase, "-READ-singleton.fastq.gz"]))
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
            try:
                shutil.rmtree(os.path.join(pth, 'split-adapter-quality-trimmed'))
            except:
                pass


class FileOptions():
    def __init__(self):
        self.tworeads = False
        self.read1 = ''
        self.read2 = ''

    def get_complex_arguments(self, conf, section='params'):
        self.tworeads = conf.getboolean('params', 'separate reads')
        self.read1 = conf.get('params', 'read1')
        self.read2 = conf.get('params', 'read2')

    def get_arguments(self, conf, section='params'):
        self.tworeads = conf.getboolean('params', 'separate reads')
        self.read1 = conf.get('params', 'read1')


def main():
    message()
    args = get_args()
    conf = ConfigParser.ConfigParser()
    conf.optionxform = str
    conf.read(args.conf)
    nproc = multiprocessing.cpu_count()
    options = FileOptions()
    if nproc >= 2 and args.cores >= 2:
        pool = multiprocessing.Pool(args.cores)
    else:
        pool = None
    if args.remap:
        names = conf.items('remap')
    else:
        names = conf.items('map')
    if args.complex and conf.has_section('params'):
        options.get_complex_arguments(conf)
    elif conf.has_section('params'):
        options.get_arguments(conf)
    create_new_dir(args.output, None)
    sample_map = get_tag_names_from_sample_file(args.input, names, args.remap, options)
    newpths = make_dirs_and_rename_files(args.input, args.output, sample_map, args.rename, args.copy, options)
    if args.only_cleanup:
        cleanup_intermediate_files(newpths, args.interleave)
        sys.exit()
    if args.adapter:
        if options.tworeads:
            np = [[i, conf, 'split-adapter-trimmed'] for i in newpths]
        else:
            np = [[i, conf, 'adapter-trimmed'] for i in newpths]
        trim_adapter_sequences(pool, np)
    if args.quality:
        if args.pe:
            if not options.tworeads:
                split_reads(pool, newpths)
            trim_low_qual_reads(pool, newpths)
    if args.pe and args.interleave:
        interleave_reads(pool, newpths)
    if args.cleanup:
        cleanup_intermediate_files(newpths, args.interleave)
    print ""

if __name__ == '__main__':
    main()
