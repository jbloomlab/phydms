.. phydms documentation master file, created by
   sphinx-quickstart on Tue Mar 31 11:35:01 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for ``phydms``
==================================

`phydms`_ is a software package for **phy**\logenetic analyses using **d**\eep **m**\utational **s**\canning data to inform the substitution models.

`phydms`_ can use :ref:`ExpCM` to describe the evolution of protein-coding genes for phylogenetic inference and the detection of biologically interesting selection. 

`phydms`_ is written in `Python`_ by `Jesse Bloom`_. The heavy lifting of the phylogenetic likelihood calculations is performed using the `Bio++`_ libraries (thanks to `Julien Dutheil`_ and `Laurent Gueguen`_ for making these wonderful libraries available). 

:ref:`installation` of `phydms`_ will install the command-line executable ``phydms``, which performs the phylogenetic analyses. It also installs the auxillary programs ``phydms_comprehensive`` and ``phydms_plotselection``, which facilitate phylogenetic model comparisons and visualization of site-specific selection.

The `phydms source code`_ is freely available under a `GPLv3`_ license (the necessary `Bio++`_ libraries are packaged into the `phydms source code`_, and are themselves under a `GPL-compatible`_ `CeCILL`_ license). However, it is recommended that you install `phydms`_ from `PyPI`_ using `pip`_ as described in the :ref:`installation` instructions.

If you use `phydms`_, please cite the references in :ref:`acknowledgments`.

Contents
------------

.. toctree::
   :maxdepth: 1

   installation
   quick_example
   phydms_prog
   phydms_comprehensive_prog
   phydms_plotselection_prog
   ExpCM
   api
   acknowledgments
   indices_and_tables

.. include:: weblinks.txt
