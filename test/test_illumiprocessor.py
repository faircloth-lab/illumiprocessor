#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 31 January 2014 12:36 PST (-0800)
"""


import os
import glob
import hashlib

#import pdb


class TestGetTruHtReads:
    def test_enough_reads(self, fake_truht_reads):
        assert len(fake_truht_reads) == 2

    def test_correct_file_names(self, fake_truht_reads):
        expected_r1 = set([
            'fake-truht_S1_L001_R1_001.fastq.gz',
            'fake-truht_S2_L001_R1_001.fastq.gz'
        ])
        expected_r2 = set([
            'fake-truht_S1_L001_R2_001.fastq.gz',
            'fake-truht_S2_L001_R2_001.fastq.gz'
        ])
        observed_r1 = []
        observed_r2 = []
        for read in fake_truht_reads:
            observed_r1.extend([os.path.basename(r) for r in read.r1])
            observed_r2.extend([os.path.basename(r) for r in read.r2])
        assert set(observed_r1) == expected_r1
        assert set(observed_r2) == expected_r2


class TestReadTrimmingResults:
    def check_read_files(self, fake_truht_args, hashes, dir):
        for file in glob.glob(os.path.join(
                fake_truht_args.output,
                dir,
                "split-adapter-quality-trimmed",
                "*.fastq.gz"
        )):
            md5 = hashlib.md5(open(file, 'rb').read()).hexdigest()
            assert md5 == hashes[os.path.basename(file)]

    def test_read_trimming(self, fake_truht_args):
        from illumiprocessor.main import main
        main(fake_truht_args)
        test_truht1 = {
            'adapters.fasta': '5868b3e17fc058bd54cb2f4d8445e289',
            'fake-truht1-READ-singleton.fastq.gz': '9c01c7ea1c5c87be1e1df4eb9ed91372',
            'fake-truht1-READ1.fastq.gz': 'dadf96e04afc519f0e375981b369b925',
            'fake-truht1-READ2.fastq.gz': '1ae68e11567670040c2d48dc5d96d9f7',
        }
        self.check_read_files(fake_truht_args, test_truht1, "test_truht1")
        test_truht2 = {
            'adapters.fasta': '1d5a8f37a1bf503743faeaf0f4263394',
            'fake-truht2-READ-singleton.fastq.gz': 'a02633175b2df700a8c8dc1a48ab67cc',
            'fake-truht2-READ1.fastq.gz': 'e513950c3d7766d04457e67405004ac4',
            'fake-truht2-READ2.fastq.gz': '1e1080a9c3bffb6e308dee1e46c5f2aa'
        }
        self.check_read_files(fake_truht_args, test_truht2, "test_truht2")
