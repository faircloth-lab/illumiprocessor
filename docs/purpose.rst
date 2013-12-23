Purpose
=======

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

These names correspond to a general file structure for each sample that looks
like::

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
downstream processing (like `phyluce`_).

This is why I wrote `illumiprocessor`_ - it processess many input files of
`Illumina`_ data, in parallel, to:

1. rename hundreds of fastq files
2. create a sample-specific, ``adapters.fasta`` file for adapter trimming
3. trim adapter contamination from input fastq files
4. trim low quality bases from input fastq files

`illumiprocessor`_ is a wrapper script around a software package written in JAVA
named `trimmomatic`_ that runs `trimmomatic`_ against many `Illumina`_ fastq
files in parallel.

In our hands, `trimmomatic`_ is the best adapter and quality trimmer we have
used, and it is developed and maintained by `Björn Usadel's group
<http://www.usadellab.org>`_. `trimmomatic`_ outperforms a number of other read
trimmers, it is reasonably fast, and it offers a lots of nice trimming options.

.. include:: links.rst
