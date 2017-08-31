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

  * *ExpCM* with a gamma-distributed :math:`\beta` (if using the ``--gammabeta`` flag).

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
A file with the suffix ``modelcomparison.md`` will be created that summarizes the model comparison.
For each model, it reports the :math:`\Delta\rm{AIC}`, the optimized log likelihood, and the values of key parameters.
This file is in Markdown format::

    | Model                   | deltaAIC | LogLikelihood | nParams | ParamValues                                   |
    |-------------------------|----------|---------------|---------|-----------------------------------------------|
    | ExpCM_NP_prefs          | 0.00     | -3389.38      | 6       | beta=2.99, kappa=6.31, omega=0.78             |
    | averaged_ExpCM_NP_prefs | 2586.44  | -4682.60      | 6       | beta=0.28, kappa=6.51, omega=0.12             |
    | YNGKP_M5                | 2599.70  | -4683.23      | 12      | alpha_omega=0.30, beta_omega=2.41, kappa=5.84 |
    | YNGKP_M0                | 2679.50  | -4724.13      | 11      | kappa=5.79, omega=0.11                        |


``phydms`` output for each model
++++++++++++++++++++++++++++++++++
For each individual model, there will also be all of the expected ``phydms`` output files as described in :ref:`phydms_prog`. These files will begin with the prefix specified by ``outprefix``, which will be followed by the name of the model.

.. include:: weblinks.txt
