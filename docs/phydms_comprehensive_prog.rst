.. _phydms_comprehensive_prog:

=======================================
the ``phydms_comprehensive`` program
=======================================

.. contents::

Overview
------------
``phydms_comprehensive`` is a program that simplifies usage of ``phydms`` for standard analyses. Essentially, ``phydms_comprehensive`` runs ``phydms`` for several different models to enable model comparisons and identify selection.

In its simplest usage, you simply provide ``phydms_comprehensive`` with an alignment and one or more files giving site-specific amino-acid preferences.
The program then first uses ``RAxML`` to infer a tree under the *GTRCAT* model.
Alternatively, you can specify the tree using ``--tree`` flag.
For each set of site-specific amino-acid preferences, ``phydms_comprehensive`` optimizes the tree with the following models:

  * *ExpCM*

  * *ExpCM* with a gamma-distributed :math:`\omega` (if using the ``--gammaomega`` flag).

  * *YNGKP_M0*

  * *YNGKP_M5*

It also runs each of these *ExpCM* with averaged preferences as a control.
Finally, it creates summaries that enable comparison among the models.

You could get all of this output by simply running ``phydms`` repeatedly, but using ``phydms_comprehensive`` automates this process for you.

See below for information on `Command-line usage`_ and `Output files`_.

Command-line usage
--------------------
.. argparse::
   :module: parsearguments
   :func: PhyDMSComprehensiveParser
   :prog: phydms_comprehensive


   outprefix
    This prefix can be a directory name (e.g. ``my_directory/``) if you want to create a new directory. See `Output files`_ for a description of the created files.

   alignment
    The alignment must meet the same specifications described in the ``phydms`` documentation for the argument of the same name (see :ref:`phydms_prog`).

   prefsfiles
    Provide the name of one or more files giving site-specific amino-acid preferences. These files should meet the same specifications described in the ``phydms`` documentation for the *prefsfile* that should accompany an *ExpCM* (see :ref:`phydms_prog`).

   \-\-raxml
     By default, ``phydms_comprehensive`` uses ``RAxML`` to infer a tree topology under the *Jukes Cantor* model and the command ``raxml``. Use this argument to specify a path to ``RAxML`` other than ``raxml``.  The tree inferred by ``RAxML`` can be found in the same directory as the other output files.

   \-\-tree
     If you want to instead fix the tree to some existing topology, use this argument and provide the name of a file giving a valid tree in Newick format. You cannot specify both ``--tree`` and ``--raxml``.

Output files
--------------
Running ``phydms_comprehensive`` will create the following output files, all with the prefix specified by ``outprefix``.

Log file
++++++++++++
A file with the suffix ``.log`` will be created that summarizes the overall progress of ``phydms_comprehensive``. If ``outprefix`` is just a directory name, this file will be called ``log.log``.

Model comparison files
+++++++++++++++++++++++++
Two model comparison files will be created, one with the suffix ``modelcomparison.csv`` and one with the suffix ``modelcomparison.md`` to summarize the model fitting.
The first summary file includes all of the parameters and the final log likelihood for each model.
The second file reports a subset of the parameters, the final log likelihood, and the :math:`\Delta` AIC in reference to the best-fit model.


``phydms`` output for each model
++++++++++++++++++++++++++++++++++
For each individual model, there will also be all of the expected ``phydms`` output files as described in :ref:`phydms_prog`. These files will begin with the prefix specified by ``outprefix``, which will be followed by a description of the model. For instance, if the command is::

    phydms_comprehensive my_directory/ alignment.fasta prefs1.txt prefs2.txt

then we expect output files with the following prefixes:

    * ``my_directory/ExpCM_prefs1_*`` : the *ExpCM* using the preferences in ``prefs1.txt``.

    * ``my_directory/averaged_ExpCM_prefs1_*`` : the *ExpCM* using the preferences in ``prefs1.txt`` averaged across sites.

    * ``my_directory/ExpCM_prefs2_*`` : the *ExpCM* using the preferences in ``prefs2.txt``.

    * ``my_directory/averaged_ExpCM_prefs2_*`` : the *ExpCM* using the preferences in ``prefs2.txt`` averaged across sites.

    * ``my_directory/YNGKP_M0_*`` : the *YNGKP_M0* model.

    * ``my_directory/YNGKP_M8_*`` : the *YNGKP_M8* model.


.. include:: weblinks.txt
