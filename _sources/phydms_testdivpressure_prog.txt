.. _phydms_testdivpressure_prog:

=======================================
the ``phydms_testdivpressure`` program
=======================================

.. contents::

Overview
------------
``phydms_testdivpressure`` is a program that simplifies usage of ``phydms`` to compare different models of diversifying pressure. 

In its simplest usage, you provide ``phydms_testdivpressure`` with an alignment, a file giving site-specific amino-acid preferences, and a file giving site-specific diversifying pressures. 
The program first uses ``phydms`` to infer a tree under the *YNGKP_M0* model (starting from a crude neighbor-joining tree) and then optimizes the tree under the *ExpCM* (without any diversifying pressures specified).
The program then fits an *ExpCM* augmented by each set of diversifying pressures to the optimized tree.
Finally, it creates a summary that enables comparison among the models.

For an additional control, ``phydms_testdivpressure`` can randomize the diversifying pressures and fit an *ExpCM* using these randomized diversifying pressures. 

If this program hangs up on the initial inference of the tree with *YNGKP_M0*, then infer a tree using another tree-building program such as `codonPhyML`_ and pass that tree via the ``--treetopology`` option.

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
    The format should be same as that described in the ``phydms`` documentation (see :ref:`phydms_prog`).

   divpressure
    The format should be the same as that described in the ``phydms`` documentation (see :ref:`phydms_prog`).

   \-\-treetopology
    By default, ``phydms_testdivpressure`` uses ``phydms`` to infer a tree topology under the *YNGKP_M0* model starting from a crude neighbor-joining tree. If you want to instead fix the tree to some existing topology, use this argument and provide the name of a file giving a valid tree in Newick format.

   \-\-randomizations
    This number specifies the number of times the diversifying pressures are randomized. An adequate number of simulations will give an approximate p-value.
    
   \-\-optimizebrlen
    By default, ``phydms_testdivpressure`` fixes the tree topology optimized under an *ExpCM* without diversifying pressure. Then for each analysis with diversifying pressure, it only re-scales all branch lengths by a single parameter (this is equivalent to passing the ``--fixbrlen`` and ``--addrateparameter`` to ``phydms``). If instead you want to optimize all the branch lengths for each diversifing pressure analysis, use this argument.

   \-\-randomseed
    This number specifies the first random seed used to randomized the diversifying pressure and the random seed number increments by one for each subsequent randomization.

Output files
--------------
Running ``phydms_testdivpressure`` will create the following output files, all with the prefix specified by ``outprefix``.

Log file
++++++++++++
A file with the suffix ``.log`` will be created that summarizes the overall progress of ``phydms_testdivpressure``. If ``outprefix`` is just a directory name, this file will be called ``log.log``.

Model comparison file
+++++++++++++++++++++++++
A file with the suffix ``_modelcomparison.csv`` (or just the name ``modelcomparison.txt`` if ``outprefix`` is just a directory name) will be created that summarizes the model fitting. Here is an example of the file that would be created if ``phydms_testdivpressure`` was passed a single ``divpressure`` with the name ``divpressure1.txt`` and ``--randomizations`` was set to one.


.. csv-table:: modelcomparison.csv
   :header: "Name","DiversifyingPressureName","DiversifyingPressureType","LogLikelihood","omega","omega2","kappa","beta","phiA","phiC","phiG"
   :widths: 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10

   "noDiversifyingpressure","none","none",-7306.12,0.709701,"N/A",4.069475,2.045563,0.388787,0.186856,0.2175507
   "divpressure1","divpressure1","true",-7295.2,0.576956,1.027609,4.097212,2.031828,0.387518,0.187517,0.217087
   "divpressure1_random_0","divpressure1","random",-7305.14,0.737622,-0.241151,4.073150,2.044787,0.388257,0.187045,0.2176979

``phydms`` output for each model
++++++++++++++++++++++++++++++++++
For each individual model, there will also be all of the expected ``phydms`` output files as described in :ref:`phydms_prog`. These files will begin with the prefix specified by ``outprefix``, which will be followed by a description of the model. For instance, if the command is::

    phydms_testdivpressure my_directory/ alignment.fasta prefs.txt divpressure1.txt divpressure2.txt --randomizations 1
    
then we expect output files with the following prefixes:

    * ``my_directory/divpressure1*`` : the *ExpCM* using the diversifying pressures in ``divpressure1.txt``.

    * ``my_directory/divpressure1_random_0*`` : the *ExpCM* using a randomized version of the diversifying pressures in ``divpressure1.txt``.

    * ``my_directory/divpressure2*`` : the *ExpCM* using the diversifying pressures in ``divpressure2.txt``.

    * ``my_directory/divpressure2_random_0*`` : the *ExpCM* using a randomized version of the diversifying pressures in ``divpressure2.txt``.



.. include:: weblinks.txt
   
