*********************
Using illumiprocessor
*********************

Using `illumiprocessor`_ is reasonably simple after the depdencies are
installed.  It requires two steps:

1. Create a configuration file telling the program what to do
2. Run the program against the configuration file

Creating a configuration file
=============================

The configuration file is structured using the standard python configuration
format where sections are denoted using brackets (``[section]``), and key-value
pairs are placed under each section given paramters (keys) and their values
(values).  The `illumiprocessor`_ configuration file has four sections:

1. The adapter structure section
2. The sequence tag section
3. The map of sequence tags to samples
4. The map of original sample names to new sample names

The entire file
---------------

The entire file structure looks like the following, which is explained in
detail, below::

    [adapters]
    i7:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
    i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

    [tag sequences]
    BFIDT-030:ATGAGGC
    BFIDT-003:AATACTT

    [tag map]
    F09-44_ATGAGGC:BFIDT-030
    F09-96_AATACTT:BFIDT-003

    [names]
    F09-44_ATGAGGC:F09-44
    F09-96_AATACTT:F09-96

[adapters]
----------

The adapter structure section is headed by::

    [adapters]

and the heading is followed by two sequences that denote the adapter
structures::

    [adapters]
    i7:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
    i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

The adapters are labelled using a **required** naming scheme - the first is
prefixed with `i7`, which, if the insert DNA strand is oriented 5' to 3', refers
to the "right" side adapter and the second is prefixed with `i5`, which, if the
insert DNA strand is oriented the same way, refers to the "left" side adapter.

We use an asterisk to denote where, in the adapter structure, the  sample-
specific index sequence is added.  **This will be filled in automatically for
each sample and sample-specific adapter**.

Below are several, indexed adapter structures commonly used with `Illumina`_
sequencing.


Illumina TruSeq v3 (single-indexed)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `Illumina`_ TruSeq v3 adapters are largely used with "older" PE library
preparations for multiplexed samples.  The `i7` adapter contains the sequence
tag (AKA "barcode" or "index"), and the system is a "single-indexed" system,
meaning that only one index is used to identify samples::

    [adapters]
    i7:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
    i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

Illumina TruSeq HT (double-indexed)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Illumina`_ TruSeq HT adapters are "double-indexed" and both the `i7` and `i5`
adapters contain the sequence tag (AKA "barcode" or "index").  Because the
system uses two indexes, both indexes are used in concert to identify samples,
and the config file needs to have two asterisks indicating where these indexes
get inserted::

    [adapters]
    i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
    i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT

Illumina TruSeq LT (single-indexed)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can convert he `Illumina`_ TruSeq HT system to what is called the TruSeq
LT system by using only one of the two indexes **on the i7 adapter**.  This
makes the TruSeq LT system functionally equivalent to the older TruSeq v3
system, except that the adapter sequences are different::

    [adapters]
    i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
    i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGTGTAGATCTCGGTGGTCGCCGTATCATT

Illumina Nextera (XT)
^^^^^^^^^^^^^^^^^^^^^

The `Illumina`_ Nextera system is another way of preparing libraries that uses
enzymatic shearing by transposase to (1) shear the DNA and (2) integrate
adapters to the DNA for subsequent sequencing.  The benefits of this system are
speed.  The `Illumina`_ Nextera system uses double-indexing, similar to the
`Illumina`_ TruSeq HT system, except that the adapter sequences are different.
The config file needs to have two asterisks indicating where these indexes get
inserted::

    [adapters]
    i7:CTGTCTCTTATACACATCTCCGAGCCCACGAGAC*ATCTCGTATGCCGTCTTCTGCTTG
    i5:CTGTCTCTTATACACATCTGACGCTGCCGACGA*GTGTAGATCTCGGTGGTCGCCGTATCATT

[tag sequences]
---------------

The tag sequence section maps a "tag name" (or set of tag names) onto the actual
tag sequences that correspond to each "tag name".  In this way, we're able to
refer to the "tag name", which is often clearer than using the tag sequence
(which is just a combination of A, C, G, & T).

Single-indexed libraries
^^^^^^^^^^^^^^^^^^^^^^^^

The format for single-indexed libraries looks like the following::

    [tag sequences]
    BFIDT-030:ATGAGGC
    BFIDT-003:AATACTT

This means that there are two libraries we're processing, one with the tag
``BFIDT-030`` and another with the tag ``BFIDT-003``.  The sequence on the
other side of the colon from each tag name will be inserted into the ``i7``
adapter at the asterisk.

The tag sequences should be input in 5' to 3' orientation.


Dual-indexed libraries
^^^^^^^^^^^^^^^^^^^^^^

The format for double-indexed libraries is very similar to the above, except
that the section generally contains more tags **and the tag names must include
the** ``i5`` **and** ``i7`` **designations**::

    [tag sequences]
    i7-N701:GCTACGCT
    i7-N702:GGACTCCT
    i5-N501:TAGATCGC
    i5-N502:CTCTCTAT

So, here, we've denoted two different sorts of tags that have names prepended with
``i7`` and ``i5``.  This means that the ``i7`` sequences are those prepended
with ``i7`` and the ``i5`` sequences are those prepended with ``i5``.

The tag sequences should be input in 5' to 3' orientation.

[tag map]
---------

The ``[tag map]`` section is where we denote which ``fastq.gz`` files output by
the sequencer contain which tags.  There are two slightly different formats for
this section depending on whether libraries are single-indexed or dual-indexed.

Single-indexed libraries
^^^^^^^^^^^^^^^^^^^^^^^^

The format for single-indexed libraries looks like the following::

    [tag map]
    morelia-viridis1_GCCTTCA:BFIDT-030
    cnemldophorus-sexlineatus1_GGTACGC:BFIDT-003

This means that the fastq.file whose name contains (note the wildcard)::

    morelia-viridis1_GCCTTCA_L005_R*_001.fastq.gz

likely contains adapter contamination, and the index used to identify those
samples is from BFIDT-030.  In other words, here, we're mapping the tag
BFIDT-030, whose sequence we denoted before, onto the sample whose fastq files
are::

    morelia-viridis1_GCCTTCA_L005_R1_001.fastq.gz
    morelia-viridis1_GCCTTCA_L005_R2_001.fastq.gz

.. note:: You should only input the fastq filename up to the end of the
   sequence tag - the remainder of the file name is filled in using a
   relatively standard wildcard structure.  If this wildcard structure does not
   fit your samples, you can run ``illumiprocessor`` using the ``--r1-pattern``
   and ``--r2-pattern`` arguments.

Dual-indexed libraries
^^^^^^^^^^^^^^^^^^^^^^

The format for dual-indexed libraries looks like the following::

    [tag map]
    morelia-viridis1_GCCTTCA:i7-N701,i5-N501
    cnemidophorus-sexlineatus1_GGTACGC:i7-N702,i5-N502

This means that the fastq.file whose name contains (note the wildcard)::

    morelia-viridis1_GCCTTCA_L005_R*_001.fastq.gz

likely contains adapter contamination, and the indexes used to identify those
samples is the combination of i7-N701 **and** i5-N501.  Here, we're mapping both
of the tags i7-N701 and i5-N501, whose sequence we denoted before, onto the
sample whose fastq files are::

    morelia-viridis1_GCCTTCA_L005_R1_001.fastq.gz
    morelia-viridis1_GCCTTCA_L005_R2_001.fastq.gz

.. note:: You should only input the fastq filename up to the end of the
   sequence tag - the remainder of the file name is filled in using a
   relatively standard wildcard structure.  If this wildcard structure does not
   fit your samples, you can run ``illumiprocessor`` using the ``--r1-pattern``
   and ``--r2-pattern`` arguments.

[names]
-------

The names section remaps `Illumina`_-specific file names onto something that's
genreally more pleasing for the end-user.  For instance, we can place the
following into the ``[names]`` section::

    morelia-viridis1_GCCTTCA:morelia-viridis1
    cnemldophorus-sexlineatus1_GGTACGC:cnemidophorus-sexlineatus1

which takes the files originally named::

    morelia-viridis1_GCCTTCA_L005_R1_001.fastq.gz
    morelia-viridis1_GCCTTCA_L005_R2_001.fastq.gz

and renames them, after trimming adapter contamination and low-quality bases
to (see the :ref:`Output` section below for more info)::

    morelia-viridis1-READ1.fastq.gz
    morelia-viridis1-READ2.fastq.gz

Running the program
===================

Once you have setup the configuration file, the program is ready to run.  You
run the program using the following::

    illumiprocessor \
        --input <path-to-directory-of-read-files-to-clean> \
        --output <path-to-directory-of-cleaned-reads-to-output> \
        --config <path-to-config-file> \
        --cores 12

.. program:: illumiprocessor

.. cmdoption:: --input (required)

   The PATH to the input data (a folder of fastq.gz files).

.. cmdoption:: --output (required)

    The PATH to where you want to store the output data.

.. cmdoption:: --config (required)

   The PATH to the config files.

.. cmdoption:: --cores (optional; default = 1)

    The number of compute cores to run simultaneously.

.. cmdoption:: --trimmomatic (optional; default = ~/anaconda/bin/trimmomatic.jar)

    The PATH to the `trimmomatic`_ jar file.

.. cmdoption:: --min-len (optional; default = 40)

    The minimum length of trimmed sequences to output.

.. cmdoption:: --no-merge (optional; default = False)

    Do not merge singleton output files.

.. cmdoption:: --r1-pattern (optional; default = None)

    A regular expression pattern for matching R1 reads.

.. cmdoption:: --r2-pattern (optional; default = None)

    A regular expression pattern for matching R2 reads.

.. cmdoption:: --se (optional; default = False)

    Trim single-end (SE) reads.

.. cmdoption:: --phred (optional; default = PHRED33)

    The quality scoring system used for the read data (PHRED33 or PHRED64).

.. cmdoption:: --log-path (optional; default = ./)

    A path to a directory in which to store the log file(s) output.

.. cmdoption:: --verbosity (optional; default = INFO)

    The verbosity level to use for log file output.


.. _Output:

The output directory structure and files
========================================

After running the program, the resulting directory structure at the requested
output path will look like::

    output-folder/
        morelia-viridis1/
            adapters.fasta
            raw-reads/
                [symlink to R1]
                [symlink to R2]
            split-adapter-quality-trimmed/
                morelia-viridis1-READ1.fastq.gz
                morelia-viridis1-READ2.fastq.gz
                morelia-viridis1-READ-singleton.fastq.gz
            stats/
                morelia-viridis1-name-adapter-contam.txt
        cnemidophorus-sexlineatus1/
            adapters.fasta
            raw-reads/
                [symlink to R1]
                [symlink to R2]
            split-adapter-quality-trimmed/
                cnemidophorus-sexlineatus1-READ1.fastq.gz
                cnemidophorus-sexlineatus1-READ2.fastq.gz
                cnemidophorus-sexlineatus1-READ-singleton.fastq.gz
            stats/
                cnemidophorus-sexlineatus1-name-adapter-contam.txt


The "singleton" file
--------------------

The ``singleton`` file in the ``split-adapter-quality-trimmed`` directory
contains all of the paired reads that lost their "mate" (the other member of
the pair) due to trimming.  `trimmomatic`_ outputs these reads separately, but
`illumiprocessor`_ combines them together in a single file.

.. include:: links.rst
