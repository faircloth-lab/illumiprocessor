Installation
=============

JAVA
-----

Illumiprocessor uses `trimmomatic`_, which is a JAVA program.  As a result, you
need to install JAVA for your platform.  This is beyond the scope of this
document, but here are some pointers.

osx
^^^

It is best to install JAVA 1.7 on osx, as that version also supports GATK and
several other useful tools.  To install JAVA 1.7, download and install the
package here: http://www.java.com/en/download/manual.jsp

centos linux
^^^^^^^^^^^^^
::

    su -c "yum install java-1.7.0-openjdk"

ubuntu linux
^^^^^^^^^^^^
::

    sudo apt-get install openjdk-7-jre

conda
------

We use the `anaconda`_ distribution with the `conda`_ package manger to manage
almost all of our package for `python`_.  There are several reasons for this,
one being that we can manage lots of **types** of packages with conda.  Another
is that `conda`_ manages dependencies of packages **very, very well**.

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
locations that will work in a number of situations.

standard
---------

Of course, you can install `illumiprocessor`_ in the standard way.  As above
you need to ensure that you have:

1. installed JAVA
2. installed `trimmomatic`_

Once those are completed, you can download the `source`_, then run::

    python setup.py install

The ``illumiprocessor`` binary will be where your `python` distribution stores
binaries, and the libraries will be in the `site-packages` directories for your
$PYTHONPATH.

.. _illumiprocessor: https://github.com/faircloth-lab/illumiprocessor
.. _trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. _anaconda: https://store.continuum.io/cshop/anaconda/
.. _conda: http://docs.continuum.io/conda/
.. _scythe: https://github.com/vsbuffalo/scythe
.. _sickle: https://github.com/najoshi/sickle
.. _documentation: http://illumiprocessor.readthedocs.org/
.. _python: http://python.org/
.. _source: https://github.com/faircloth-lab/illumiprocessor/releases

