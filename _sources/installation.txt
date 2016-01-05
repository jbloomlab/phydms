.. _installation:

================================
Installation
================================

.. contents::

Minimal requirements
----------------------
Installing `phydms`_ requires `Python`_ version 2.7, the `Python`_ package management system `pip`_, and a ``g++``/``c++`` compiler. These instructions assume you are using a reasonably modern Linux operating system, although it may also be possible to install `phydms`_ on Mac OS X (this has not been carefully checked). If you have all of this, any other dependencies should be automatically installed if you follow the `Quick installation`_ instructions.

Quick installation
---------------------
If you have the `Minimal requirements`_, you should be able to install `phydms`_ by simply typing::

    pip install phydms --user

or::

    sudo pip install phydms

depending on whether you are `Installing locally versus globally`_.

Note that installation will take a little while (maybe 5-15 minutes) as all of the `Bio++`_ source code included with `phydms`_ must be compiled.

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
You can also install from the `phydms source code`_ on GitHub. Again, if you are simply using `phydms`_ (and not developing it), you are strongly suggested to **not** do a source installation and instead use `pip`_.

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
`phydms`_ requires `Python`_ 2.7 and a ``g++``/``c++`` compiler. 

In addition it requires some external `Python`_ packages. The up-to-date exact requirements are listed in the ``setup.py`` file in the main directory of the `phydms source code`_. If you are installing with `pip`_, these external packages will automatically be installed. If you are installing from source, you will need to install these packages yourself. As of the writing of this documentation, these packages are:

    * `Biopython`_ is used for sequence and tree manipulations; you need at least version 1.65

    * `cython`_ is used to build the ``c++`` extensions. You need at least version 0.21.

    * `dms_tools`_ is used for some file input / output operations, as well as to create the plots created by ``phydms_plotselection``. `dms_tools`_ itself has dependencies that include `weblogo`_ and `matplotlib`_, but these should automatically be installed if you are using `pip`_ to `phydms`_ (or `dms_tools`_). You need `dms_tools`_ version at least 1.1.5.

    * `scipy`_ is used for optimization and some statistical / mathematical operations. You need at least version 0.16.

    * `matplotlib`_ is used for plotting. You need at least version 1.3.

What about ``Bio++``?
------------------------
Most of the actual calculations performed by `phydms`_ use the `Bio++`_ libraries, which are themselves written in ``c++``. Manual installation of `Bio++`_ and subsequent linking to the installation of `phydms`_ is somewhat complicated, so the relevant portions of the `Bio++`_ source code are simply included in the `phydms source code`_ (in ``phydms/phydmslib/Bpp/`` you will find ``bpp-core``, ``bpp-seq``, and ``bpp-phyl``). When you install `phydms`_ using `pip`_ (or using ``python setup.py install``), `cython`_ compiles this `Bio++`_ source code in a form that is usable by `phydms`_. The requirement for compiling `Bio++`_ is why it takes a little while to install `phydms`_.

If you are developing `phydms`_, you might get tired of constantly re-compiling all of `Bio++`_. Or if you are testing new developments in `Bio++`_, you might want to link `phydms`_ to a separate installation of `Bio++`_.

If you are installing from source code, the ``setup.py`` script includes an option that should allow you to instead dynamically link to pre-compiled `Bio++`_ libraries. To do this, you need to first install the correct versions of ``bpp-core``, ``bpp-seq``, and ``bpp-phyl`` to the expected (standard) location for global or local installation. Then run::

    python setup.py install --user --dynamically-link-bpp

or::

    python setup.py install --dynamically-link-bpp

depending on whether you are `Installing locally versus globally`_. The ``--dynamically-link-bpp`` option tells ``setup.py`` to instead link to these pre-installed `Bio++`_ libraries. If this isn't working (e.g. you don't have the `Bio++`_ libraries in the location expected by ``setup.py``), you will have to delve into the ``setup.py`` file yourself to figure out what is going wrong.

.. include:: weblinks.txt
