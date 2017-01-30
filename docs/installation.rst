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
If you have the `Minimal requirements`_, you should be able to install `phydms`_ by simply typing::

    pip install phydms --user

or::

    sudo pip install phydms

depending on whether you are `Installing locally versus globally`_.

Installing locally versus globally
---------------------------------------
You need to figure out whether you want to install `phydms`_ globally or locally. 

Install globally
~~~~~~~~~~~~~~~~~~~~

If you are the administrator (super-user) for your computer, then you can install programs globally. In that case, you will typically need to run the install commands prefaced with ``sudo``, and will be asked for your password. 

Install locally
~~~~~~~~~~~~~~~~~~~~
If you are using a shared computing cluster, then you are probably not the administrator. In that case, you need to install `phydms`_ locally in a directory that you own. The `Python convention for local installation`_ is to install in ``~/.local/`` by adding the option ``--user`` to installation commands.

For locally installed programs to be accessible, you need to add ``~/.local/bin/`` to the ``PATH`` variable, and ``~/.local/lib/`` to the ``PYTHONPATH`` variable. If you are using the `bash shell`_, you would do this by adding the following lines to your ``~/.bashrc`` file::

    PATH=$HOME/.local/bin/:$PATH
    export PYTHONPATH=$HOME/.local/lib/python2.7:$PATH

You then want to make sure that your ``~/.bash_profile`` simply sources your ``~/.bashrc`` file as `described here <http://www.joshstaiger.org/archives/2005/07/bash_profile_vs.html>`_ by making ``~/.bash_profile`` consist of the following::

    if [ -f ~/.bashrc ]; then
        source ~/.bashrc
    fi

More details on installing with ``pip``
----------------------------------------------------
The advantage of using `pip`_ is that it will automatically take care of installing the `Other required software`_. Here are some more details of installing with `pip`_ if the `Quick installation`_ fails.

First, make sure ``pip`` is installed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check if you have `pip`_ installed. You can do this by typing::

    pip -h

If you get the help message for `pip`_, then `pip`_ is installed and you can move to the next step. If you instead get an error message such as ``-bash: pip: command not found`` then you need to install `pip`_.

If you have ``easy_install``, then you can install `pip`_ globally with::

    sudo easy_install pip

or locally with::

    easy_install pip --user

If those commands also fail (i.e. you don't have ``easy_install`` either), then install `pip`_ by following the `instructions here <https://pip.pypa.io/en/latest/installing.html>`_.

Next, use ``pip`` to install ``phydms``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once `pip`_ is installed, you can do a local installation with::

    pip install phydms --user

or a global installation with::

    sudo pip install phydms

Upgrading your ``phydms`` version with ``pip``
--------------------------------------------------
If you have previously installed `phydms`_ but are not sure that you have the latest version, you can upgrade using `pip`_. To do this for a local installation, use::

    pip install phydms --user --upgrade

or globally::

    sudo pip install phydms --upgrade

Install ``phydms`` from source code
-----------------------------------------------------------------------------------
You can also install from the `phydms source code`_ on GitHub. Again, if you are simply using `phydms`_ (and not developing it), you are suggested to instead simply use `pip`_.

To install from source, first make sure that you have the `Other required software`_.

Then download or clone the `phydms source code`_ from GitHub.
After unpacking the main ``phydms`` directory that contains this source, install locally with::

    cd phydms
    python setup.py install --user

or globally with::

    cd phydms
    sudo python setup.py install 

If you want to install from source but want `pip`_ to process the dependencies, you can install locally with::

    cd phydms
    pip install -e . --user

or globally with::

    cd phydms
    pip install -e .

Other required software
------------------------------------------
`phydms`_ requires some external `Python`_ packages. The up-to-date exact requirements are listed under ``install_requires`` in the ``setup.py`` file in the main directory of the `phydms source code`_. If you are installing with `pip`_, these external packages will automatically be installed. If you are installing from source, you will need to install these packages yourself. 


.. include:: weblinks.txt
