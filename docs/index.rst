.. illumiprocessor documentation master file, created by
   sphinx-quickstart on Sun Dec 22 12:14:07 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

illumiprocessor: parallel adapter and quality trimming
======================================================

Release v\ |version|. (:ref:`Changelog`)

`illumiprocessor`_ is a tool to batch process illumina sequencing reads using
the excellent `trimmomatic`_ package. The program takes as input the location
of your demultiplexed illumina reads, a configuration file that is formatted
in Microsoft Windows INI file format (key:value pairs, see the example file),
and the output directory in which you want to store the trimmed reads.

`illumiprocessor`_ will trim adapter contamination and low quality bases from
SE and PE illumina reads, and the program is also capable of dealing with
double-indexed reads . The current version of `illumiprocessor`_ uses
`trimmomatic`_ instead of `scythe`_ and `sickle`_ (used in v1.x) because we
have found the performance of trimmomatic to be better, particularly when
dealing with double-indexed illumina reads.  However, you may find that
running `scythe`_ after trimming with `illumiprocessor`_ may ensure that every
bit of potential adapter contamination is removed.

`illumiprocessor`_ is suited to parallel processing in which each set of
illumina reads are processed on a separate (physical) compute core.
`illumiprocessor`_ **assumes that all fastq files input to the program
represent individuals samples (i.e., that you have merged mulitple files for
each read from the same sample by combining fastq.gz files)**.


Guide
=====

.. toctree::
   :maxdepth: 2

   purpose
   install
   usage
   disclaimers


Project info
============

.. toctree::
   :maxdepth: 1

   citing
   license
   changelog
   authors
   contributors


.. include:: links.rst
