.. _installation:

================================
Installation
================================

.. contents::

Minimal requirements
----------------------
`phydms`_ is written in `Python`_. You can run it either under ``python2`` (if you have version 2.7 or higher) or ``python3`` (if you have version 3.4 or higher).

Straightforward installation requires the `Python`_ package management system `pip`_ and a ``C`` compiler such a ``gcc`` (there are some ``cython`` extensions). 

`phydms`_ has been tested on relatively recent versions of Linux and Mac OS X.

Quick installation
---------------------
If you have the `Minimal requirements`_, you can install `phydms`_ by simply typing::

    pip install phydms --user

This assumes that the local installation directory is in your current path. See below for more details.

Where to install `phydms`_ with ``pip``
-----------------------------------------
You need to figure out where you want to install `phydms`_.
Using ``sudo`` to install globally `is not recommended for Python packages in general <http://stackoverflow.com/questions/21055859/what-are-the-risks-of-running-sudo-pip/21056000#21056000>`_.

We suggest that you use one of the following two options:

Install via the ``--user`` option to ``pip``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The `Python convention for local installation`_ is to install locally using the ``--user`` option to ``pip``. This option installs into ``~/.local/`` for Linux, or into ``$HOME/.local/lib/Python/x.y/`` (where ``x.y`` is the `Python`_ version, such as ``3.5``) for Mac OS X.

For locally installed programs to be accessible, you need to add ``~/.local/bin/`` to the ``PATH`` variable, and ``~/.local/lib/`` to the ``PYTHONPATH`` variable. If you are using the `bash shell`_, you would do this by adding the following lines to your ``~/.bashrc`` file::

    PATH=$HOME/.local/bin/:$PATH
    export PYTHONPATH=$HOME/.local/lib/python3.5/:$PATH

You then want to make sure that your ``~/.bash_profile`` simply sources your ``~/.bashrc`` file as `described here <http://www.joshstaiger.org/archives/2005/07/bash_profile_vs.html>`_ by making ``~/.bash_profile`` consist of the following::

    if [ -f ~/.bashrc ]; then
        source ~/.bashrc
    fi

Once the paths are set up as described above, simply install with::

    pip install --user

If this fails, check if you have `pip`_ installed. You can do this by typing::

    pip -h

If you instead get an error message such as ``-bash: pip: command not found`` then you need to install `pip`_.
If you have ``easy_install``, then you can install `pip`_ globally with::

    sudo easy_install pip

or locally with::

    easy_install pip --user

If those commands also fail (i.e. you don't have ``easy_install`` either), then install `pip`_ by following the `instructions here <https://pip.pypa.io/en/latest/installing.html>`_.

Install into a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The other good approach is to use ``pip`` to install into a virtual environment. 
The use of virtual environments `is described here <http://stackoverflow.com/questions/21055859/what-are-the-risks-of-running-sudo-pip/21056000#21056000>`_. Once you have created a virtual environment as in these instructions, you can just install with ``pip`` dropping the ``--user`` command.

Upgrading your ``phydms`` version with ``pip``
--------------------------------------------------
If you have previously installed `phydms`_ but are not sure that you have the latest version, you can upgrade using `pip`_::

    pip install phydms --upgrade --user


Install ``phydms`` from source code
-----------------------------------------------------------------------------------
First clone the `phydms source code`_ from GitHub::

    git clone https://github.com/jbloomlab/phydms


Then install locally with::

    cd phydms
    pip install -e . --user

Other required software
------------------------------------------
`phydms`_ requires some external `Python`_ packages. The up-to-date exact requirements are listed under ``install_requires`` in the ``setup.py`` file in the main directory of the `phydms source code`_. If you are installing with `pip`_, these external packages will automatically be installed. If you are installing from source, you may need to install these packages yourself. 

Older versions (e.g., version 1.0)
-----------------------------------
The original version of `phydms`_ (version ``1.*``) used in `Bloom, Biology Direct, 12:1`_ had similar functionality, but an entirely different implementation that utilized `Bio++`_ for the likelihood calculations.
The new version (version ``2.*``) has similar functionality, but is a completely new implementation and so may give slightly different results.

The last stable version of the old implementation is ``1.3.0``. 
The source code for that version is `available here <https://github.com/jbloomlab/phydms/tree/1.3.0>`_.
The older implementations are also still on `PyPI`_.

.. include:: weblinks.txt
