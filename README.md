# illumiprocessor [![Build Status](https://travis-ci.org/faircloth-lab/illumiprocessor.png?branch=master)](https://travis-ci.org/faircloth-lab/illumiprocessor)

[illumiprocessor][1] is a tool to batch process illumina sequencing reads using
the excellent [trimmomatic][2] package. The program takes a configuration file
that is formatted in Microsoft Windows INI file format (key:value pairs, see the
example file).

[illumiprocessor][1] will trim adapter contamination from SE and PE illumina
reads and is capable of dealing with double-indexed reads and read trimming
(example to come shortly). The current version of [illumiprocessor][1] uses
[trimmomatic][2] instead of [scythe][3] and [sickle][4] (used in v1.x) because
we have found the performance of trimmomatic to be better, particularly when
dealing with double-indexed illumina reads.  However, you may find that running
[scythe][3] after trimming with [illumiprocessor][1] or [trimmomatic][2] ensures
that every bit of potential adapter contamination is removed.

[illumiprocessor][1] is suited to parallel processing in which each set of
illumina reads are processed on a separate (physical) compute core.
**[illumiprocessor][1] assumes that all fastq files input to the program
represent individuals samples (i.e., that you have merged mulitple files for
each read from the same sample by combining fastq.gz files)**.

## Citing Illumiprocessor

If you use illumiprocessor in your work, you can cite the software as follows:

    Faircloth, BC. 2013. illumiprocessor: a trimmomatic wrapper for parallel
    adapter and quality trimming. http://dx.doi.org/10.6079/J9ILL.

Please be sure also to cite [trimmomatic][2]:

    Lohse M, Bolger AM, Nagel A, Fernie AR, Lunn JE, Stitt M, Usadel B. 2012.
    RobiNA: a user-friendly, integrated software solution for RNA-Seq-based
    transcriptomics. Nucleic Acids Res. 40(Web Server issue):W622-7.
    doi:10.1093/nar/gks540

    Del Fabbro C, Scalabrin S, Morgante M and Giorgi FM. 2013. An Extensive
    Evaluation of Read Trimming Effects on Illumina NGS Data Analysis. PLoS
    ONE 8(12): e85024. doi:10.1371/journal.pone.0085024


## installation

Illumiprocessor uses [trimmomatic][2], which is a JAVA program, **so you need 
to install JAVA for your platform**.

### conda

If you are using [anaconda][5] or the [conda][6] package manager, you can
automatically install everything you need by editing `~/.condarc` to add the
`faircloth-lab` repository, so that the file looks like:

    # channel locations. These override conda defaults, i.e., conda will
    # search *only* the channels listed here, in the order given. Use "default"
    # to automatically include all default channels.

    channels:
      - defaults
      - http://conda.binstar.org/faircloth-lab

Then run:

    conda install illumiprocessor

This will install [trimmomatic][2] and [illumiprocessor][1].

### standard

Ensure that you have installed JAVA.  Install [trimmomatic][2].  Once those are
completed, download the source, then:

    python setup.py install

## quick-start

To run [illumiprocessor][1], you setup a config file (`<path-to-config-file>`)
like:

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

Then you run illumiprocessor against the config file using:

    illumiprocessor \
        --input <path-to-directory-of-read-files-to-clean> \
        --output <path-to-directory-of-cleaned-reads-to-output> \
        --config <path-to-config-file> \
        --cores 12

This will output a directory containing reads organised using the following
structure:

    sample1-name/
        adapters.fasta
        raw-reads/
            [symlink to R1]
            [symlink to R2]
        split-adapter-quality-trimmed/
            sample1-name-READ1.fastq.gz
            sample1-name-READ2.fastq.gz
            sample1-name-READ-singleton.fastq.gz
        stats/
            sample1-name-adapter-contam.txt
    sample2-name/
        ...
    sample3-name

# More information

For more information and a more complete description of all of these steps,
please see the [documentation][7].

[1]: https://github.com/faircloth-lab/illumiprocessor
[2]: http://www.usadellab.org/cms/?page=trimmomatic
[3]: https://github.com/vsbuffalo/scythe
[4]: https://github.com/najoshi/sickle
[5]: https://store.continuum.io/cshop/anaconda/
[6]: http://docs.continuum.io/conda/
[7]: http://illumiprocessor.readthedocs.org/
