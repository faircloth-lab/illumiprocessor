The purpose of illumiprocessor
==============================

`Illumina`_ sequencing offers incredible capacity.  Often, biologists wish to
split this capacity among many samples run in multiplex.  When you do this, and
you demultiplex the samples that you run together, you generally get files
output from `Casava`_ that look something like::

    acutotyphlops-kunuaensis_CCGGTGG_L005_R1_001.fastq.gz
    acutotyphlops-kunuaensis_CCGGTGG_L005_R2_001.fastq.gz

    amphisbaena-fuliginosa_GCTCCTC_L005_R1_001.fastq.gz
    amphisbaena-fuliginosa_GCTCCTC_L005_R2_001.fastq.gz

    boiga-irregularis_GAGAGTA_L005_R1_001.fastq.gz
    boiga-irregularis_GAGAGTA_L005_R2_001.fastq.gz

These names correspond to a general structure that is::

    <name>_<sequence-tag>_<lane>_<read>_<read-file>.fastq.gz

where the ``<name>`` is equivalent to the species binomial name used above, the
string of letters is equivalent to ``<sequence-tag>``, the ``<lane>`` the
samples were run in was Lane 5 (L005), each sample has a ``<read>`` Read1 (R1)
and Read2 (R2) file because this was paired-end (PE) sequencing, and each file
has only a single ``<read-file>`` (001).

When you multiplex many samples into the same lane, you can end up with many
many output files (e.g. 100).  But, we still need to trim adapter contamination
and low quality bases from these files.  Adapter contamination occurs when the
insert DNA sequence is shorter than the read length, and we sequence part of the
adapter during the sequencing process.  Low quality bases occur as a function
of read length, chemistry, the sequencing platform, etc.  Removal of adapter
contamination and low quality bases is essential to proper downstream analysis
and processing.

The problem is that with so many files, removing adapter contamination and low
quality bases is problematic across all of the files.  Doing so may require
tedious processing by hand (which is bad) or shell scripting, which many people
are not comfortable with.  When we trim adapters and remove contamination, we
may also want to do things such as bulk-renaming of the sequence files and
placing the resulting trimmed data in a static directory structure for
downstream processing with a package like `phyluce`_.

This is why I originally `illumiprocessor`_ - it processess input files of
`Illumina`_ data, in parallel to:

1. rename hundreds of fastq files
2. trim adapter contamination from input fastq files
3. trim low quality bases from input fastq files

`illumiprocessor`_ is a parallel wrapper script around a software package
written in JAVA named `trimmomatic`_.  `trimmomatic`_ is the best adapter and
quality trimmer we have used and it is developed and maintained by `Björn
Usadel's group <http://www.usadellab.org>`_.  In my hands, `trimmomatic`_
outperforms a number of other read trimmers, it is reasonably fast, and it
offers a lots of nice options.

.. include:: links.rst
