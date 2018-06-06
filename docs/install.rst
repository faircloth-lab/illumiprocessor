.. include:: links.rst

************
Installation
************

The instructions below assume you are either:

1. running osx
2. running linux (x86_64)

We do not support any other platforms, although you should be able to install
and run `illumiprocessor`_ and its dependencies on various flavors of windows.

trimmomatic & illumiprocessor
=============================

You can install the illumiprocessor_ dependencies in one of two ways:  (1)
using `conda`_ and (2) manually.  We **strongly suggest** that you use
`conda`_.  There are several reasons for this, one being that we can manage
lots of **types** of packages with conda. Another is that `conda`_ manages
dependencies of packages **very, very well**.

.. attention:: Manual installation of illumiprocessor_ is not supported.

.. attention:: illumiprocessor is already installed as part of phyluce_. You
    do not need to install anything else after installing phyluce_.

First, you need to install anaconda_ or miniconda_ **with Python 2.7**.  Whether
you choose miniconda_ or anaconda_ is up to you, your needs, how much disk
space you have, and if you are on a fast/slow connection.

.. attention:: You can easily install anaconda_ or miniconda_ in your ``$HOME``,
    although you should be aware that this setup can sometimes cause problems in
    cluster-computing situations.

.. tip:: Do I want anaconda_ or miniconda_?
    :class: admonition tip

    The major difference between the two python distributions is that anaconda_
    comes with many, many packages pre-installed, while miniconda_ comes with
    almost zero packages pre-installed.  As such, the beginning anaconda_
    distribution is roughly 200-500 MB in size while the beginning miniconda_
    distribution is 15-30 MB in size.

    **We suggest that you install miniconda.**

.. tip:: What version of miniconda_ or anaconda_ do I need?
    :class: admonition tip

    Right now, illumiprocessor_ **only runs with Python 2.7**.  This means that you need
    to install a version of miniconda_ or anaconda_ that uses Python 2.7.  The
    easiest way to do this is to choose carefully when you download a
    particular distribution for your OS (be sure to choose the Python 2.7
    version).

miniconda
^^^^^^^^^

Follow the instructions here for your platform:
https://conda.io/docs/user-guide/install/index.html

.. note:: Once you have installed either Miniconda or Anaconda, we will refer
    to the install as `conda` throughout the remainder of this documentation.

anaconda
^^^^^^^^

Follow the instructions here for your platform:
http://docs.continuum.io/anaconda/install.html


.. note:: Once you have installed either Miniconda or Anaconda, we will refer
    to the install as `conda` throughout the remainder of this documentation.


Checking your `$PATH`
^^^^^^^^^^^^^^^^^^^^^

Regardless of whether you install miniconda_ or anaconda_, you need to check
that you've installed the package correctly.  To ensure that the correct
location for anaconda_ or miniconda_ are added to your $PATH (this occurs
automatically on the $BASH shell), run the following::

    $ python -V

The output should look similar to (`x` will be replaced by a version)::

    Python 2.7.x :: Anaconda x.x.x (x86_64)

Notice that the output shows we're using the `Anaconda x.x.x` version of
Python_. If you do not see the expected output (or something similar), then you
likely need to edit your $PATH variable to add anaconda_ or miniconda_.

The easiest way to edit your path, if needed is to open ``~/.bashrc`` with a
text editor (if you are using ZSH, this will be ``~/.zshrc``) and add, as the
last line::

    export PATH=$HOME/path/to/conda/bin:$PATH

where ``$HOME/path/to/conda/bin`` is the location of anaconda/miniconda on your
system (usually ``$HOME/anaconda/bin`` or ``$HOME/miniconda/bin``).

.. warning:: If you have previously set your ``$PYTHONPATH`` elsewhere in your
   configuration, it may cause problems with your anaconda_ or miniconda_
   installation of illumiprocessor_.  The solution is to remove the offending library
   (-ies) from your ``$PYTHONPATH``.


Add the necessary bioconda repositories to conda
------------------------------------------------

You need to add the location of the bioconda_ repositories to your conda_
installation.  To do that, you can follow the instructions `at the bioconda
site <https://bioconda.github.io/#set-up-channels>`_ or you can simply edit
your ``~/.condarc`` file to look like:

.. code-block:: bash

    channels:
      - defaults
      - conda-forge
      - bioconda

Once you do this, you have access to all of the packages installed at
bioconda_ and conda-forge_.  The order of this file is important - conda_ will
first search in it's default repositories for package, then it will check
conda-forge, finally it will check bioconda.

How to install illumiprocessor
-----------------------------

You now have two options for installing illumiprocessor_.  You can install
illumiprocessor_ in what is known as a `conda environment
<https://conda.io/docs/user-guide/tasks /manage-environments.html>`_, which
lets you keep code for different applications separated into different
environments.  **We suggest this route**.

You can also install all of the illumiprocessor_ code and dependencies in
your default conda_ environment.

Install illumiprocessor in it's own conda environment
-----------------------------------------------------

We can install everything that we need for illumiprocessor_ in it's own environment by running:

.. code-block:: bash

    conda create --name illumiprocessor illumiprocessor

This will create an environment named ``illumiprocessor``, then download and install
everything you need to run illumiprocessor_ into this `illumiprocessor` conda environment. To
use this illumiprocessor environment, you **must** run:

.. code-block:: bash

    source activate illumiprocessor

To stop using this illumiprocessor environment, you **must** run:

.. code-block:: bash

    source deactivate

Install illumiprocessor in the default conda environment
-------------------------------------------------------

We can simply install everything that we need in our default conda_
environment, as well.  In some ways, this is easier, but it could be viewed as
a less-ideal option in terms of repeatability and separability of functions.
To install illumiprocessor_ in the default environment, after making sure that you
have miniconda_ or anaconda_ in your $PATH, and after adding the bioconda
repositories, run:

.. code-block:: bash

    conda install illumiprocessor