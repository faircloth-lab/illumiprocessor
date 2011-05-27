
#!/usr/bin/env python
# encoding: utf-8

"""
pre_process_reads.py

Created by Brant Faircloth on 26 May 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.

DESCRIPTION
===========

Pre-process Illumina reads for assembly by renaming, removing adapter
contamination, quality trimming, removing N-bases, giving some summary
info on tags, and zipping everything up.

Assumes that Bin_005.zip has been extracted into a dir of the same name e,g.:

Bin_005/s_4_sequence.txt
Bin_005/s_5_sequence.txt
Bin_005/s_6_sequence.txt
Bin_005/s_7_sequence.txt

Based on contents of SampleDirectories.csv, it will rename the input files
w/ their respective "names" as given in SampleDirectories.csv.

USAGE
=====

python process_reads.py --sample-map SamplesDirectories.csv Bin_005

"""

import os
import sys
import glob
import zlib
import shutil
import argparse
import multiprocessing
from collections import Counter
from tools.sequence import fastq
from tools.sequence import transform

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Parse fastq files and drop reads containing Ns.')
    parser.add_argument('input', help='the input directory containing the reads')
    parser.add_argument('--sample-map', required = True, dest = 'sample_map')
    return parser.parse_args()

def get_tag_names_from_sample_file(sample_directories):
    sample_map = {}
    for k,v in enumerate(open(sample_directories, 'rU')):
        if k == 0:
            pass
        else:
            values = v.strip().split(',')
            lane, tag, species, bin = values[1], values[4], values[5], values[-1]
            sample_map['_'.join((bin, lane))] = [tag, species]
    return sample_map

def make_dirs_and_rename_files(sample_map, tld):
    [os.system('mkdir -p {0}'.format(os.path.join(tld, i))) for i in ['raw', 'adapt-trim', 'qual-trim', 'n-less', 'stats']]
    for name in glob.glob(os.path.join(tld,'*.txt')):
        temp_bin, temp_name = os.path.split(name)
        bin_num = temp_bin.split('_')[1]
        lane_num = temp_name.split('_')[1]
        real_name = sample_map['_'.join((bin_num, lane_num))][1].replace(' ', '-') + '.fastq'
        #pdb.set_trace()
        shutil.move(name, os.path.join(temp_bin, 'raw', real_name))

def scythe_runner(reads):
    tld, filename = reads
    inpt = os.path.join(tld, "raw", filename)
    outp = os.path.join(tld, "adapt-trim", filename)
    log = os.path.join(tld, "adapt-trim", "{0}.log".format(filename))
    os.system('scythe -a adapters.fasta -o {0} {1} >> {2}'.format(outp, inpt, log))
    sys.stdout.write(".")
    sys.stdout.flush()
    return

def trim_adapter_sequences(pool, tld, reads = 'raw'):
    sys.stdout.write("\nTrimming adapter contamination")
    fastqs = [[tld, os.path.basename(f)] for f in glob.glob(os.path.join(tld, reads,'*.fastq'))]
    pool.map(scythe_runner, fastqs)
    return

def sickle_runner(reads):
    tld, filename = reads
    inpt = os.path.join(tld, "adapt-trim", filename)
    outp = os.path.join(tld, "qual-trim", filename)
    log = os.path.join(tld, "qual-trim", "{0}.log".format(filename))
    os.system("sickle se -f {0} -t illumina -o {1} -q 20 -l 40 > {2} 2>&1".format(inpt, outp, log))
    sys.stdout.write(".")
    sys.stdout.flush()
    return

def trim_low_qual_reads(pool, tld, reads = 'adapt-trim'):
    sys.stdout.write("\nTrimming low quality reads")
    fastqs = [[tld, os.path.basename(f)] for f in glob.glob(os.path.join(tld, reads,'*.fastq'))]
    pool.map(sickle_runner, fastqs)
    return

def drop_n_reads_runner(reads):
    tld, filename = reads
    inpt = os.path.join(tld, "qual-trim", filename)
    outp = os.path.join(tld, "n-less", filename)
    in_fastq = fastq.FastqReader(inpt)
    out_fastq = fastq.FastqWriter(outp)
    for read in in_fastq:
        if 'N' not in read.sequence:
            out_fastq.write(read)
        else:
            pass
    out_fastq.close()
    sys.stdout.write(".")
    sys.stdout.flush()
    return

def drop_n_reads(pool, tld, reads = 'qual-trim'):
    sys.stdout.write("\nDropping reads containing n-bases")
    fastqs = [[tld, os.path.basename(f)] for f in glob.glob(os.path.join(tld, reads,'*.fastq'))]
    pool.map(drop_n_reads_runner, fastqs)
    return

def get_sequence_tags_runner(reads):
    tld, filename = reads
    inpt = os.path.join(tld, "raw", filename)
    outp = open(os.path.join(tld, "stats", filename), 'w')
    cnt = Counter()
    for read in fastq.FastqReader(inpt):
        cnt[read.identifier.split('#')[-1].split('/')[0]] += 1
    for tag in cnt.keys():
        outp.write("{0},{1},{2}\n".format(tag, cnt[tag], transform.DNA_reverse_complement(tag)))
    outp.close()
    sys.stdout.write(".")
    sys.stdout.flush()
    return

def get_sequence_tags(pool, tld, reads = 'raw'):
    sys.stdout.write("\nGetting sequence tag counts")
    fastqs = [[tld, os.path.basename(f)] for f in glob.glob(os.path.join(tld, reads,'*.fastq'))]
    pool.map(get_sequence_tags_runner, fastqs)
    return

def zip_shit_up_runner(reads):
    tld, filename = reads
    inpt = os.path.join(tld, "n-less", filename)
    z = open(os.path.join(tld, "n-less", filename + '.gz'), 'wb')
    contents = open(inpt, 'rb').read()
    z_contents = zlib.compress(contents)
    z.write(z_contents)
    z.close()
    sys.stdout.write(".")
    sys.stdout.flush()
    return

def zip_shit_up(pool, tld, reads = 'n-less'):
    sys.stdout.write("\nZipping shit up")
    fastqs = [[tld, os.path.basename(f)] for f in glob.glob(os.path.join(tld, reads,'*.fastq'))]
    pool.map(zip_shit_up_runner, fastqs)
    return

def main():
    args = get_args()
    # careful here - zlib calls eat RAM
    pool = multiprocessing.Pool(4)
    sample_map = get_tag_names_from_sample_file(args.sample_map)
    make_dirs_and_rename_files(sample_map, args.input)
    get_sequence_tags(pool, args.input)
    trim_adapter_sequences(pool, args.input)
    trim_low_qual_reads(pool, args.input)
    drop_n_reads(pool, args.input)
    zip_shit_up(pool, args.input)
    print ""
    
if __name__ == '__main__':
    main()
