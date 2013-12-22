# illumiprocessor

[illumiprocessor][1] is a tool to batch process illumina sequencing reads using the
excellent [trimmomatic][2] package. The program takes a configuration file that is
formatted in Microsoft Windows INI file format (key:value pairs, see the example file).

[illumiprocessor][1] will trim adapter contamination from SE and PE illumina reads and is
capable of dealing with double-indexed reads and read trimming (example to come shortly).
The current version of [illumiprocessor][1] uses [trimmomatic][2] instead of [scythe][3]
and [sickle][4] (used in v1.x) because we have found the performance of trimmomatic to be
better, particularly when dealing with double-indexed illumina reads.  However, you may find
that running [scythe][3] after trimming with [illumiprocessor][1] or [trimmomatic][2] ensures
that every bit of potential adapter contamination is removed.

[illumiprocessor][1] is suited to parallel processing in which each set of illumina reads
are processed on a separate (physical) compute core.  **[illumiprocessor][1] assumes that all
fastq files input to the program represent individuals samples (i.e., that you have
merged mulitple files for each read from the same sample by combining fastq.gz files)**.

To run [illumiprocessor][1], you setup a config file (`<path-to-config-file>`) like:

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

This will output a directory containing reads organised using the following structure:

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

[1]: https://github.com/faircloth-lab/illumiprocessor
[2]: http://www.usadellab.org/cms/?page=trimmomatic
[3]: https://github.com/vsbuffalo/scythe
[4]: https://github.com/najoshi/sickle
