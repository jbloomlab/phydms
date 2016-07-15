.. _phydms_comprehensive_prog:

=======================================
the ``phydms_comprehensive`` program
=======================================

.. contents::

Overview
------------
``phydms_comprehensive`` is a program that simplifies usage of ``phydms`` for standard analyses. Essentially, ``phydms_comprehensive`` runs ``phydms`` for several different models to enable model comparisons and identify selection.

In its simplest usage, you simply provide ``phydms_comprehensive`` with an alignment and one or more files giving site-specific amino-acid preferences. The program then first uses ``phydms`` to infer a tree under the *YNGKP_M0* model (starting from a crude neighbor-joining tree), and then for each set of site-specific amino-acid preferences optimizes this tree (with fixed topology) for an *ExpCM* and a control *ExpCM* with the preferences averaged across sites (the ``phydms`` option ``--avgprefs``). The program also optimizes the tree using the *YNGKP_M8* model. For all models, it detects site-specific selection. It also creates a simple summary that enables comparison among the models.

If this program hangs up on the initial inference of the tree with *YNGKP_M0*, then infer a tree using another tree-building program such as `codonPhyML`_ and then pass that via the ``--treetopology`` option.

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

   \-\-treetopology
    By default, ``phydms_comprehensive`` uses ``phydms`` to infer a tree topology under the *YNGKP_M0* model starting from a crude neighbor-joining tree. If you want to instead fix the tree to some existing topology, use this argument and provide the name of a file giving a valid tree in Newick format.

   \-\-no-omegabysite 
    By default, ``phydms_comprehensive`` infers a site-specific :math:`\omega_r` for the tree obtained using each *YNGKP* model and each *ExpCM*. Use this option if you do **not** want to do this.

   \-\-no-stringencybysite 
    By default, ``phydms_comprehensive`` infers a site-specific :math:`\beta_r` for each *ExpCM*. Use this option if you do **not** want to do this.

   \-\-no-diffprefsbysite 
    By default, ``phydms_comprehensive`` infers site-specific :math:`\Delta\pi_{r,a}` for each *ExpCM*. Use this option if you do **not** want to do this.

   \-\-diffprefconc
    See the documentation for ``phydms`` for more details.

   \-\-no-avgprefs
    As described in the documentation for ``phydms`` (see :ref:`phydms_prog`), there are two sensible controls when using an *ExpCM* model: averaging the preferences across sites or randomizing them across sites. By default, ``phydms_comprehensive`` performs the averaging control (this is the ``--avgprefs`` option to ``phydms``). Use this option if you do **not** want to include this averaging control.

   \-\-yngkp
    In addition to the *YNGKP_M0* model (which is always used to optimize the initial tree), the program will also optimize any additional *YNGKP* model variants listed here. See the documentation for :ref:`phydms_prog` for more information about these model variants.

    If you don't want to include any additional *YNGKP* models beyond *M0*, use ``--yngkp ""``.

   \-\-randprefs
    As described in the documentation for ``phydms`` (see :ref:`phydms_prog`), there are two sensible controls when using an *ExpCM* model: averaging the preferences across sites or randomizing them across sites. By default, ``phydms_comprehensive`` performs the averaging control (this is the ``--avgprefs`` option to ``phydms``) but not the randomization control (this is the ``--randprefs`` option to ``phydms``). Use this option if you also want to include the ``--randprefs`` control.

   \-\-avgrandcontrol
    As described in the documentation for ``phydms`` (see :ref:`phydms_prog`), there are two sensible controls for an *ExpCM* model: randomizing preferences among sites and averaging them across sites. The ``--no-avgprefs`` and ``--randprefs`` options specify whether to perform these controls for the *ExpCM* informed by each of the potentially numerous preferences specified in files listed in ``prefsfiles``. If you have a lot of preference files, you might only want to perform these controls for one of them. In that case, list that one file (which should also be listed in ``prefsfiles``) here. The only control will then be average and random preferences for this one file.

   \-\-use_existing
    Perhaps some of the output that ``phydms_comprehensive`` will create already exists, such as from a previous run. By default, ``phydsm_comprehensive`` deletes any existing output and re-creates everything from scratch. If you don't want to do that, then use this option. But be careful: ``phydms_comprehensive`` only looks for existing output by file name; if those files were actually generated using a different set of input data or program options, ``phydms_comprehensive`` has no way of determining this.

Output files
--------------
Running ``phydms_comprehensive`` will create the following output files, all with the prefix specified by ``outprefix``.

Log file
++++++++++++
A file with the suffix ``.log`` will be created that summarizes the overall progress of ``phydms_comprehensive``. If ``outprefix`` is just a directory name, this file will be called ``log.log``.

Model comparison file
+++++++++++++++++++++++++
A file with the suffix ``_modelcomparison.txt`` (or just the name ``modelcomparison.txt`` if ``outprefix`` is just a directory name) will be created that summarizes the model fitting. Here is an example of the file that would be created if ``phydms_comprehensive`` is passed a single ``prefsfiles`` with the name ``prefs.txt``::

    ======================== ====== ============== ===========================================================================================
    model                    AIC    log likelihood number parameters (optimized + empirical): optimized values                                
    ======================== ====== ============== ===========================================================================================
    ExpCM_prefs              0.0    -4415.2        6 (6 + 0): beta = 2.92, omega = 0.77, kappa = 6.19, phiA = 0.37, phiC = 0.20, phiG = 0.22  
    YNGKP_M8                 2482.1 -5648.3        14 (5 + 9): pomegas = 0.08, omegas = 1.02, betap = 0.00, betaq = 4.34, kappa = 6.40
    averaged_ExpCM_prefs     2600.9 -5715.7        6 (6 + 0): beta = 0.32, omega = 0.12, kappa = 6.51, phiA = 0.35, phiC = 0.18, phiG = 0.25  
    YNGKP_M0                 2619.1 -5719.8        11 (2 + 9): omega = 0.11, kappa = 6.24                                                     
    ======================== ====== ============== ===========================================================================================

In this table, ``ExpCM_prefs`` is the *ExpCM* created using the preferences in ``prefs.txt``, and ``averaged_ExpCM_prefs`` is the *ExpCM* created after averaging the preferences in ``prefs.txt`` across sites (the ``--avgprefs`` option to ``phydms``. 

For the *YNGKP_M8* model, ``pomegas`` is the weight assigned to the positively selected (:math:`\omega \gt 1`) category, and ``omegas`` is the value of :math:`\omega` for this category. The ``betaq`` and ``betap`` parameters are the shape of the beta distribution for :math:`\omega < 1`.

The models are ranked by `AIC`_ in terms of difference from the best fitting model. The third column shows the number of parameters. The *YNGKP* models have 9 empirical parameters corresponding to the codon frequencies estimated using the *F3X4* option.

Note that the table is in `reStructuredText`_ format, and so can be converted to attractive HTML using ``rst2html``.

If ``--dateseqs`` is being used, then a file giving the estimated dates to most recent common ancestor (MRCA) is also generated with the suffix ``_mrca_dates.txt`` (or just the name ``mrca_dates.txt`` if ``outprefix`` is just a directory name). Here is an example of such a file::

    ======================= =========
    model                   MRCA date
    ======================= =========
    ExpCM NP prefs          1911.14
    YNGKP M0                1911.87
    averaged ExpCM NP prefs 1911.95
    YNGKP M8                1912.05
    ======================= =========

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
   
