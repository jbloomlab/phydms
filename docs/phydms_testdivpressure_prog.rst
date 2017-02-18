.. _phydms_testdivpressure_prog:

=======================================
the ``phydms_testdivpressure`` program
=======================================

.. contents::

Overview
------------
``phydms_testdivpressure`` is a program that simplifies usage of ``phydms`` to compare different models of diversifying pressure.

In its simplest usage, you simply provide ``phydms_testdivpressure`` with an alignment, a file giving site-specific amino-acid preferences, and a file giving site-specific diversifying pressures.
The program then first uses ``RAxML`` to infer a tree under the *GTRCAT* model.
Alternatively, you can specify the tree using ``--tree`` flag.
Then, the tree is optimized for a control *ExpCM* without any diversifying pressures and an *ExpCM* for each of the site-specific diversifying pressures.
For additional controls, ``phydms_testdivpressure`` will also optimize the tree with randomized diversifying pressures.
It also creates simple summary files that enable comparison among the models.

See below for information on `Command-line usage`_ and `Output files`_.

Command-line usage
--------------------
.. argparse::
   :module: parsearguments
   :func: PhyDMSTestdivpressureParser
   :prog: phydms_testdivpressure


   outprefix
    This prefix can be a directory name (e.g. ``my_directory/``) if you want to create a new directory. See `Output files`_ for a description of the created files.

   alignment
    The alignment must meet the same specifications described in the ``phydms`` documentation for the argument of the same name (see :ref:`phydms_prog`).

   prefsfile
    Provide the name a file giving site-specific amino-acid preferences. This file should meet the same specifications described in the ``phydms`` documentation for the *prefsfile* that should accompany an *ExpCM* (see :ref:`phydms_prog`).

   divpressure
    Provide the name of one or more files giving site-specific diversifying pressure values.
    The file must specify a diversifying pressure for every site in ``alignment``, using sequential 1, 2, ... numbering.
    The diversifying pressure file can be in either a comma-, tab-, or space-deliminted file with the first column giving the site and the second column giving the diversifying pressure value.
    These diversifying pressures will be scaled to fall between 0 and 1. For more information on this scaling procedure see (see :ref:`ExpCM`)

   \-\-raxml
    By default, ``phydms_comprehensive`` uses ``RAxML`` to infer a tree topology under the *Jukes Cantor* model and the command ``raxml``. Use this argument to specify a path to ``RAxML`` other than ``raxml``.  The tree inferred by ``RAxML`` can be found in the same directory as the other output files.

   \-\-tree
    If you want to instead fix the tree to some existing topology, use this argument and provide the name of a file giving a valid tree in Newick format. You cannot specify both ``--tree`` and ``--raxml``.

   \-\-randomizations
    This number specifies the number of times the diversifying pressures are randomized. An adequate number of simulations will give an approximate p-value.

Output files
--------------
Running ``phydms_testdivpressure`` will create the following output files, all with the prefix specified by ``outprefix``.

Log file
++++++++++++
A file with the suffix ``.log`` will be created that summarizes the overall progress of ``phydms_testdivpressure``. If ``outprefix`` is just a directory name, this file will be called ``log.log``.

Model comparison files
+++++++++++++++++++++++++
Two model comparison files will be created, one with the suffix ``modelcomparison.csv`` and one with the suffix ``modelcomparison.md`` to summarize the model fitting.
The first summary file includes all of the parameters of each model.
The second file summarizes the main results for the non-randomized models.

Output for individual models
++++++++++++++++++++++++++++++++++
For each individual model, there will also be all of the expected ``phydms`` output files as described in :ref:`phydms_prog`. These files will begin with the prefix specified by ``outprefix``, which will be followed by a description of the model. For instance, if the command is::

    phydms_testdivpressure my_directory/ alignment.fasta prefs.txt divpressure1.txt, divpressure2.txt --randomizations 1

then we expect output files with the following prefixes:

    * ``my_directory/ExpCM_prefs_*`` : the *ExpCM* without diversifying pressure.

    * ``my_directory/ExpCM_prefs_divpressure1_True_*`` : the *ExpCM* using the diversifying pressures in ``divpressure1.txt``.

    * ``my_directory/ExpCM_prefs_divpressure1_Random_0*`` : the *ExpCM* using a randomized version of the diversifying pressures in ``divpressure1.txt``.

    * ``my_directory/ExpCM_prefs_divpressure2_True_*`` : the *ExpCM* using the diversifying pressures in ``divpressure2.txt``.

    * ``my_directory/ExpCM_prefs_divpressure2_Random_0*`` : the *ExpCM* using a randomized version of the diversifying pressures in ``divpressure2.txt``.




.. include:: weblinks.txt
