.. phydms documentation master file, created by
   sphinx-quickstart on Tue Mar 31 11:35:01 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for ``phydms``
==================================

`phydms`_ is a software package for **phy**\logenetic analyses using **d**\eep **m**\utational **s**\canning data to inform the substitution models.

`phydms`_ can use :ref:`ExpCM` to describe the evolution of protein-coding genes for phylogenetic inference and the detection of biologically interesting selection. 

`phydms`_ is written in `Python`_ by `Jesse Bloom`_. The heavy lifting of the phylogenetic likelihood calculations is performed using the `Bio++`_ libraries (thanks to `Julien Dutheil`_ and `Laurent Gueguen`_ for making these wonderful libraries available). 

:ref:`installation` of `phydms`_ will install :ref:`phydms_prog`, which performs the phylogenetic analyses. It also installs :ref:`auxiliary_progs` that facilitate model comparisons and visualization of the results. All of these programs are installed as command-line executables.

The `phydms source code`_ is freely available under a `GPLv3`_ license (source code for `Bio++`_ and `LSD`_ is packaged into the `phydms source code`_; these codes bases are themselves under a `GPL-compatible`_ `CeCILL`_ license and a `GPLv3`_ license, respectively). 

However, you will have an easier time getting all the dependencies correct if you install `phydms`_ from `PyPI`_ using `pip`_ as described in the :ref:`installation` instructions rather than installing from the source code.

If you use `phydms`_, please cite the references in :ref:`acknowledgments`.

Contents
------------

.. toctree::
   :maxdepth: 1

   installation
   quick_example
   phydms_prog
   auxiliary_progs
   ExpCM
   api
   acknowledgments
   indices_and_tables

.. include:: weblinks.txt
