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
Then for each set of site-specific amino-acid preferences optimizes this tree (with fixed topology) for an *ExpCM* and an *ExpCM* with :math:`omega` following a :math:`Gamma` distribution.
It also runs each of these *ExpCM* models with averaged preferences as a control.  
It also creates a simple summary that enables comparison among the models.

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
Here are examples of the files that would be created if ``phydms_comprehensive`` is passed a single ``prefsfiles`` with the name ``prefs.txt``.

``modelcomparison.csv``::

  Model,parameter,value
  YNGKP_M0,kappa,5.78484
  YNGKP_M0,omega,0.11308
  YNGKP_M0,phi0A,0.340656
  YNGKP_M0,phi0C,0.153395
  YNGKP_M0,phi0G,0.308294
  YNGKP_M0,phi1A,0.32050300000000004
  YNGKP_M0,phi1C,0.206698
  YNGKP_M0,phi1G,0.236844
  YNGKP_M0,phi2A,0.331485
  YNGKP_M0,phi2C,0.202991
  YNGKP_M0,phi2G,0.243997
  YNGKP_M0,LogLikelihood,-4724.13
  ExpCM_prefs,beta,2.99308
  ExpCM_prefs,kappa,6.31356
  ExpCM_prefs,omega,0.7807689999999999
  ExpCM_prefs,phiA,0.372278
  ExpCM_prefs,phiC,0.196416
  ExpCM_prefs,phiG,0.224367
  ExpCM_prefs,LogLikelihood,-3389.38


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
