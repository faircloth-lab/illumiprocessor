************
Installation
************
`illumiprocessor`_ requires the installation of JAVA, `trimmomatic`_, and
`illumiprocessor`_.  The instructions below assume you are either:

1. running osx
2. running linux (x86_64 or osx)

We do not support any other platforms, although you should be able to install
and run `illumiprocessor`_ and its dependencies on various flavors of windows.

JAVA
====

Illumiprocessor uses `trimmomatic`_, which is a JAVA program.  As a result, you
need to install JAVA for your platform.  Getting JAVA installed is tricky across
various platforms and largely beyond the scope of this document, but here are
some very general directions for installing **JAVA 1.7**, which is best suited
for `illumiprocessor`_, along several other codebases that we use (e.g. current
GATK).

Apple OS X
-----------

To install JAVA 1.7, download and install the package from Oracle here:
http://www.java.com/en/download/manual.jsp

CentOS linux
------------

You can install the JRE with the following `yum` command::

    su -c "yum install java-1.7.0-openjdk"

Ubuntu linux
------------

You can install the JRE with the following `apt-get` command::

    sudo apt-get install openjdk-7-jre


trimmomatic & illumiprocessor
=============================

You can install the remaining dependencies in one of two ways:  (1) using
`conda`_ and (2) manually.  We **strongly suggest** that you use `conda`_, which
is part of the `anaconda`_ python distribution.  There are several reasons for
this, one being that we can manage lots of **types** of packages with conda.
Another is that `conda`_ manages dependencies of packages **very, very well**.

The use of `anaconda`_ in our lab is the default.

Installation using conda
------------------------

If you are using `anaconda`_ and/or the `conda`_ package manager, you can
automatically install everything you need by editing ``~/.condarc`` to add the
`faircloth-lab` binstar repository, so that the file looks like::

    # channel locations. These override conda defaults, i.e., conda will
    # search *only* the channels listed here, in the order given. Use "default"
    # to automatically include all default channels.

    channels:
      - defaults
      - http://conda.binstar.org/faircloth-lab

Then run::

    conda install illumiprocessor

This will install `trimmomatic`_ and `illumiprocessor`_ in standardized
locations that will work in a number of situations.  It will also test the
installations to ensure they were installed corrected.

The `trimmomatic`_ `jar` file is installed into your `$ANACONDA_HOME/jar`.


"Manual" installation
---------------------

Of course, you can install `illumiprocessor`_ in the standard (manual) way.  As
above you need to ensure that you have:

1. installed JAVA
2. installed `trimmomatic`_

Once those are completed, you can download the `illumiprocessor`_, then run::

    python setup.py install

The ``illumiprocessor`` binary will be where your `python` distribution stores
binaries, and the libraries will be in the `site-packages` directory for your
$PYTHONPATH.

.. caution:: You will need to manually specify the path to `trimmomatic`_ when
   you run `illumiprocessor`_

.. include:: links.rst
