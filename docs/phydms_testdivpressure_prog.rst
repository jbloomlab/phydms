.. _phydms_testdivpressure_prog:

=======================================
the ``phydms_testdivpressure`` program
=======================================

.. contents::

Overview
------------
``phydms_testdivpressure`` is a program that simplifies usage of ``phydms`` to compare different models of diversifying pressure. 

In its simplest usage, you simply provide ``phydms_testdivpressure`` with an alignment, a file giving site-specific amino-acid preferences, and a file giving site-specific diversifying pressures. The program then first uses ``phydms`` to infer a tree under the *YNGKP_M0* model (starting from a crude neighbor-joining tree), and then for each set of site-specific diversifying pressures optimizes this tree (with fixed topology) for an *ExpCM* and control *ExpCM* without any diversifying pressures.
For additional controls, ``phydms_testdivpressure`` will also randomize the diversifying pressures and then optimize the tree (with fixed topology) for an *ExpCM* using these randomized diversifying pressures. 
It also creates a simple summary that enables comparison among the models.

If this program hangs up on the initial inference of the tree with *YNGKP_M0*, then infer a tree using another tree-building program such as `codonPhyML`_ and then pass that via the ``--treetopology`` option.

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

   prefsfiles
    Provide the name of one or more files giving site-specific amino-acid preferences. These files should meet the same specifications described in the ``phydms`` documentation for the *prefsfile* that should accompany an *ExpCM* (see :ref:`phydms_prog`).

   \-\-treetopology
    By default, ``phydms_testdivpressure`` uses ``phydms`` to infer a tree topology under the *YNGKP_M0* model starting from a crude neighbor-joining tree. If you want to instead fix the tree to some existing topology, use this argument and provide the name of a file giving a valid tree in Newick format.

   \-\-randomizations
    This number specifies the number of times the diversifying pressures are randomized. An adequate number of simulations will give an approximate p-value.
    
Output files
--------------
Running ``phydms_testdivpressure`` will create the following output files, all with the prefix specified by ``outprefix``.

Log file
++++++++++++
A file with the suffix ``.log`` will be created that summarizes the overall progress of ``phydms_testdivpressure``. If ``outprefix`` is just a directory name, this file will be called ``log.log``.

++++++++++++++++++++++++++++++++++
For each individual model, there will also be all of the expected ``phydms`` output files as described in :ref:`phydms_prog`. These files will begin with the prefix specified by ``outprefix``, which will be followed by a description of the model. For instance, if the command is::

    phydms_testdivpressure my_directory/ alignment.fasta prefs.txt divpressure1.txt, divpressure2.txt --randomizations 1
    
then we expect output files with the following prefixes:

    * ``my_directory/divpressure1*`` : the *ExpCM* using the diversifying pressures in ``divpressure1.txt``.

    * ``my_directory/divpressure1_random_0*`` : the *ExpCM* using a randomized version of the diversifying pressures in ``divpressure1.txt``.

    * ``my_directory/divpressure2*`` : the *ExpCM* using the diversifying pressures in ``divpressure2.txt``.

    * ``my_directory/divpressure2_random_0*`` : the *ExpCM* using a randomized version of the diversifying pressures in ``divpressure2.txt``.



.. include:: weblinks.txt
   
