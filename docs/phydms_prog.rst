.. _phydms_prog:

=======================================
the ``phydms`` program
=======================================

.. contents::

Overview
------------
``phydms`` can perform phylogenetic analyses with :ref:`ExpCM` as well as with some standard non-site-specific substitution models (variants of the *YNGKP* models described in `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_). In addition, ``phydms`` can be used to detect site-specific diversifying or differential selection using the :ref:`ExpCM`.

``phydms`` is written in `Python`_, and uses the `Bio++`_ libraries to perform the core likelihood functions. See below for information on the `Command-line usage`_ and `Output files`_.

Command-line usage
--------------------
.. argparse::
   :module: parsearguments
   :func: PhyDMSParser
   :prog: phydms

   alignment
    This file should contain aligned DNA codon sequences. Stop codons are **not** allowed, except if a stop codon is the terminal character in all sequences, in which case they are automatically trimmed.

    No checking is done to remove identical sequences, so you may want to remove such sequences yourself if there a lot of them.

    The headers must be unique strings that do **not** contain any of the following: spaces, commas, colons, semicolons, parentheses, square brackets, and single or double quotes.

   tree
    This argument specifies the fixed tree topology (if **not** using ``--infertopology``) or an initial topology for further optimization. (if using ``--infertopology``). The possibilities are:

        * Specify a Newick file giving an existing tree (tip names must match sequence headers in ``alignment``). For instance, this might be a tree that you inferred using `codonPhyML`_.

        * Use *nj* to build a neighbor-joining tree from the nucleotide sequences using a crude identity model to compute distances. Trees built using this option will not be very good, so you are suggested to use this option only to get an initial tree for further optimization via ``--infertopology``.

        * Use *random* for a randomly chosen starting tree. If you use this option, you must also use ``--infertopology`` since a random tree makes no sense otherwise.

   model
    *YNGKP_M<n>* is a codon model from `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_. The *<n>* indicates the model variant. Codon frequencies are set by *F3X4* (the product of nucleotide frequency parameters at each of the three positions; see ``--fitF3X4`` for how these 9 parameters are set). The transition-transversion ratio :math:`\kappa` is optimized by maximum likelihood. The nonsynonymous-synonymous ratio :math:`\omega` is set differently depending on the variant:
     
        - *YNGKP_M0* : a single :math:`\omega` is optimized by maximum likelihood (1 free parameter)
 
        - *YNGKP_M1* : an :math:`\omega < 1` and an :math:`\omega = 1` with the weights optimized by maximum likelihood (two free parameters). This is actually the variant indicated as *M1a* in `Wong et al, Genetics, 168:1041-1051`_ rather than the original *M1* model in `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_.

        - *YNGKP_M2* : an :math:`\omega < 1`, and :math:`\omega = 1`, and an :math:`\omega > 1` with the weights optimized by maximum likelihood (4 free parameters). This is actually the variant indicated as *M2a* in `Wong et al, Genetics, 168:1041-1051`_ rather than the original *M2* model in `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_.

        - *YNGKP_M3* : three :math:`\omega` values and the weights (probabilities) of their categories optimized by maximum likelihood (5 free parameters).
 
        - *YNGKP_M7* : ``--ncats`` values of :math:`\omega` between zero and one drawn from a beta distribution with the two distribution parameters optimized by maximum likelihood (2 free parameters). 
 
        - *YNGKP_M8* : ``--ncats`` values of :math:`\omega` between zero and one drawn from a beta distribution with the two distribution parameters optimized by maximum likelihood, plus another value of :math:`\omega > 1` with the value and weight optimized by maximum likelihood (4 free parameters). 

    *ExpCM_<prefsfile>* is an :ref:`ExpCM` with amino-acid preferences taken from the file ``prefsfile``. 
    The preferences file should be in the `dms_tools`_ `preferences file format`_ for **amino acids** (any stop codon preferences if present are normalized away to zero). The preferences file must specify a preference for the amino acid encoded by every site in ``alignment``, using sequential 1, 2, ... numbering.

   outprefix
    If this prefix contains a directory name, that directory is created if it does not already exist. 

    See `Output files`_ for a description of the files created with this prefix. Any existing files with these names specified by ``outprefix`` are deleted at the start of the program's execution.

   \-\-omegabysite
    If using a YNGKP model, then the :math:`\omega_r` value is nearly analogous that obtained using the *FEL* model described by `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_. If using and *ExpCM*, then :math:`\omega_r` has the meaning described in :ref:`ExpCM`. Essentially, we fix all other model / tree parameters and then compare a model that fits a synonymous and nonsynonymous rate to each site to a null model that only fits a synonymous rate; there is evidence for :math:`\omega_r \ne 1` if fitting both nonsynonymous and synonymous rate gives sufficiently better likelihood than fitting synonymous rate alone. See also the ``--omegabysite_fixsyn`` option.

   \-\-omegabysite_fixsyn
    This option is meaningful only if you are using ``--omegabysite``. If you use this option, then we compare a model in which we fit a nonsynonymous rate to each site to a model in which we fit nothing. The synonymous rate is not fit, and so is assumed to be equal to the overall value fit for the tree. According to `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_, in some cases this can yield greater power if there is relatively limited data. However, it comes with the risk of giving spurious results if there is substantial variation in the synonymous substitution rate among sites.  

   \-\-gammarates
    Rather than a single substitution rate, the likelihood is a linear combination of those computed using four different rates drawn from a 4-category discrete gamma distribution. The shape parameter of this distribution is optimized by maixmum likelihood (the mean is fixed to one). This adds one parameter, and will increase the program's runtime.

   \-\-useLog
    This option specifies that we have `Bio++`_ use logarithms in the likelihood calculations. According to the `Bio++`_ developers, this prevents likelihoods of zero (log likelihoods of negative infinity) on large trees, but comes at the cost of slowing down the computations.

    Right now, even if ``--useLog`` is **not** used, ``phydms`` will automatically use logarithms if a specific problem is encountered while fitting site-specific differential preferences, but otherwise does not use them. Unless you have a good reason to do otherwise, keeping this setup (i.e. keeping ``--useLog`` at its default value of *False*) is probably the best choice.

   \-\-stringencybysite
    The meaning of the :math:`\beta_r` values is described in :ref:`ExpCM`.

   \-\-diffprefsbysite
    The meaning of the :math:`\Delta\pi_{r,a}` values is described in :ref:`ExpCM`.

   \-\-randprefs
    Only for *ExpCM* models. This option randomly reassigns the preferences among sites. This can be used as a control -- randomizing the preferences should make them lose their efficacy for describing evolution.

   \-\-avgprefs
    Only for *ExpCM* models. This option computes an average of each preference across sites (:math:`\pi_a = \frac{1}{L} \sum_r \pi_{r,a}` where :math:`r = 1, \ldots, L`), and then uses these average preferences for all sites. This can be used as a control, as it merges all the information in the preferences into a non-site-specific model.

   \-\-diffprefconc
    This option specifies the parameters :math:`C_1` and :math:`C_2` that determine the concentration of the prior that regularizes the differential preferences inferred by ``--diffprefsbysite``. Both these parameters must be > 0; larger values favor smaller differential preferences. See :ref:`ExpCM` for more information on the exact meaning of these parameters.

   \-\-addrateparameter
    It only makes sense to use this parameter if you have fixed all branch lengths and then wish to fit a rate that effectively scales these branch lengths.

   \-\-infertopology
    Note that sometimes the topology inference appears to fail (or at least hang up for a very long time). In that case, you may want to infer the topology with another program such as `codonPhyML`_ and then pass that topology via the ``tree`` argument rather than using ``--infertopology``.

   \-\-dateseqs
    This option specifies that we perform least-squares dating of the sequences. This may be useful if your sequences are sampled from different dates. The dating is performed using code from `LSD`_ (see `To et al, Systematic Biology, 65:82-97`_) which is embedded in `phydms`_.

    This option estimates the dates of all internal nodes and re-roots the tree. The branch lengths are adjusted in the process. Essentially, this option is equivalent to running `LSD`_ with::

        lsd -i intree -d dates -v -s seqlength -c -r as -o outtree
    
    With these options, variance weights are calculated based on the sequence length, temporal constraints are applied, and all root positions are searched with constraints on all branches.

    The file specified by this argument should give the date for each sequence in ``alignment``. The first column should be the numerical date, and the second column should give the header (exactly matching that in ``alignment``). Empty line or lines beginning with ``#`` are ignored. For instance::

        #DATE SEQUENCE
        1968.5 A/Aichi/2/1968 (H3N2) NP
        1918.0 A/Brevig Mission/1/1919 (H1N1) NP
        2006.5 A/Solomon Islands/3/2006 (H1N1) NP

   \-\-minbrlen
    Regardless of the method used to set ``tree``, all branches with lengths less than this value will be set to this value in the initial starting tree. Branches can still end up with lengths less than this after subsequent optimization of this starting tree.

   \-\-seed
    The random number seed can influence the outcome when using a non-deterministic algorithm.

   \-\-fitF3X4
    Only meaningful if ``model`` is one of the *YNGKP* models. It specifies how the 9 nucleotide frequency parameters are determined. By default, they are set to the empirically observed values in the sequence alignment. If you use ``--fitF3X4`` then these 9 parameters are fit by maximum likelihood.


   \-\-no_optimize
    This argument is intended for the case when you want to use a per-site optimization (such as ``--stringencybysite``, ``--omegabysite``, or ``--diffprefsbysite``) without re-optimizing the entire tree. It assumes that you have **previously** run ``phydms`` with the same ``outprefix`` on the same set of sequences, and have already generated the ``_modelparams.txt``, ``_tree.newick``, and ``_loglikelihood.txt`` output files. It then simply uses those files and proceeds to the per-site optimization. This option currently does **not** work for *YNGKP* models that don't use ``-fitF3X4``. When using this option, you should specify a ``treefile`` that is the one created by the previous run with ``outprefix``.

Output files
--------------
Running ``phydms`` will create the following output files, all with the prefix specified by ``outprefix``.

Log file
+++++++++++++
This is a text file with the suffix ``.log`` that records information about the program's progress.

Log likelihood file
++++++++++++++++++++++
This file has the suffix ``_loglikelihood.txt``. It simply gives the optimized log likelihood. Here is an example of a file's contents::

    log likelihood = -4415.24
                
Model parameters
++++++++++++++++++++
This file has the suffix ``_modelparams.txt``. It gives the value of all **optimized** model parameters (it does not give the value of non-optimized model parameters such as the codon frequencies under *F3X4*). Here is an example of the contents for an *ExpCM* model::

    123_Full.theta = 0.411593
    123_Full.theta1 = 0.626366
    123_Full.theta2 = 0.524156
    123_K80.kappa = 6.19349
    omega = 0.768227
    stringencyparameter = 2.92

These parameters correspond to those described in :ref:`ExpCM` as follows:

    * *stringencyparameter* is :math:`\beta`

    * *omega* is :math:`\omega`

    * *123_K80.kappa* is :math:`\kappa`

    * *123_Full.theta*, *123_Full.theta1*, and *123_Full.theta2* define the values of the nucleotide frequences :math:`\phi_A`, :math:`\phi_C`, :math:`\phi_G`, and :math:`\phi_T`. The definitions are :math:`\phi_C + \phi_G = \rm{\textit{123_Full.theta}}`, :math:`\phi_A / \left(\phi_A + \phi_T\right) = \rm{\textit{123_Full.theta1}}`, and :math:`\phi_G / \left(\phi_G + \phi_C\right) = \rm{\textit{123_Full.theta2}}`.

Here is an example of the contents for a *YNGKP_M8* model::

    kappa = 3.355859
    omegas = 3.172927
    p = 0.000100
    p0 = 0.886637
    q = 0.258658

These parameters correspond to those described for the *M8* model by `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_ (three discrete categories of :math:`\omega`), with:

    * *kappa* being :math:`\kappa`

    * *p0* being the weight assigned to the :math:`\omega < 1` beta distributed values.
   
    * *p* and *q* being the shape parameters of the beta distribution.

    * *omegas* being the :math:`\omega` value of the category (weight *1 - p0*) with :math:`\omega > 1`.
                        
Tree file
++++++++++++
This file has the suffix ``_tree.newick``, and gives the optimized tree in Newick format.

Site-specific omega file
+++++++++++++++++++++++++++
This file has the suffix ``_omegabysite.txt``, and is created only if using the ``--omegabysite`` option. This file gives the :math:`\omega_r` values optimized for each site. In the case of a *YNGKP* model, these are site-specific dN/dS ratios that should be essentially analogous to those obtained under the *FEL* model described by `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_. In the case of an *ExpCM* model, these values indicate diversifying selection for nonsynonymous amino-acid change as described in :ref:`ExpCM`.

Here is an example of the first few lines of a file::

    # Omega fit to each site after fixing tree and and all other parameters.
    # Fits compared to null model of omega = 1.000; P-values NOT corrected for multiple testing.
    #
    # site  omega   P   dLnL
    472 99.000  0.00199 4.778
    284 0.001   0.00434 4.067
    319 0.001   0.00723 3.607
    20  99.000  0.00852 3.460
    373 99.000  0.00936 3.376
    350 0.001   0.0126  3.113
    465 99.000  0.0149  2.966
    470 99.000  0.0158  2.911
    289 0.001   0.016   2.899
    396 0.001   0.024   2.546
    129 0.001   0.0253  2.501
    84  0.001   0.0279  2.418
    247 0.001   0.029   2.383
    101 8.381   0.0325  2.287

The columns should be self-explanatory, the *P*-values are for rejection of the null hypothesis that :math:`\omega_r = 1` (calculated using a :math:`\chi_1^2` test; no correction for multiple testing is included, so you need to do that yourself if necessary for your question). The sites are sorted in the file by *P*-value.

Site-specific stringency file
++++++++++++++++++++++++++++++++
This file has the suffix ``_stringencybysite.txt``, and is created only if using and *ExpCM* with ``--stringencybysite``. It gives the ratio :math:`\beta_r / \beta` described in :ref:`ExpCM`; values of this ratio greater than one indicates increased stringency for the preferences, and values less than one indicate decreased stringency for the preferences. The *P*-values are for rejection of the null hypothesis that :math:`1 = \beta_r / \beta` (calculated using a :math:`\chi_1^2` test; no correction for multiple testing is included, so you need to do that yourself if necessary for your question). The sites are sorted in the file by *P*-value. Here is an example of the first few lines of a file::

    # Stringency fit to each site after fixing tree and and all other parameters.
    # Fits compared to null model of stringency = 2.920 (the overall gene value).
    # P-values NOT corrected for multiple testing.
    # The stringency ratio is the ratio of the fitted value to the null (overall gene value).
    #
    # site  stringency_ratio    P   dLnL
    375 0.017   2.45e-05    8.903
    52  0.017   0.000249    6.709
    442 0.079   0.00048 6.097
    294 0.017   0.00399 4.145
    398 0.095   0.00858 3.454
    50  0.134   0.0103  3.291
    465 0.334   0.0123  3.132
    351 0.420   0.0125  3.121
    283 0.017   0.0126  3.110
    344 0.430   0.0287  2.392
    85  6.849   0.0338  2.253
    497 0.195   0.0341  2.244
                                        
Differential preferences file
++++++++++++++++++++++++++++++++
This file has the suffix ``_diffprefsbysite.txt`` and is created only if using *ExpCM* with ``--diffprefsbysite``. It gives the :math:`\Delta\pi_{r,a}` values for all sites and amino acids. Positive values indicate unexpectedly strong natural selection for amino acid :math:`a` at site :math:`r`, and negative values indicate unexpectedly strong selection against :math:`a` at :math:`r`.

The file is `dms_tools`_ `differential preferences file format`_. Here is an example of the first few lines::

    # POSITION WT RMS_dPI dPI_A dPI_C dPI_D dPI_E dPI_F dPI_G dPI_H dPI_I dPI_K dPI_L dPI_M dPI_N dPI_P dPI_Q dPI_R dPI_S dPI_T dPI_V dPI_W dPI_Y
    1 ? 0.0054602 -0.00109088 -0.00056453 -0.00108891 -0.00125535 -0.000617395 -0.000825462 -0.00129802 -0.00107999 -0.00249272 -0.00202186 0.0236839 -0.00145916 -0.00101602 -0.00112551 -0.00222488 -0.00140373 -0.0018838 -0.000760446 -0.000350133 -0.00112508
    2 ? 0.000905382 -0.00365815 0.00010107 0.000118787 2.0868e-05 0.000105629 9.45005e-05 0.000101899 0.000182448 0.000111788 0.000119843 2.56786e-05 0.000107759 0.000111716 0.000107354 0.000119253 0.000107919 0.000269405 0.00165618 9.42552e-05 0.000101797
    3 ? 0.0117273 -0.00227315 -1.84925e-05 0.000103205 9.00981e-05 -0.0211034 -2.02e-05 -0.0042788 4.69967e-05 6.14462e-06 -0.00705743 7.46465e-06 -0.000664213 -0.00013711 -1.53571e-05 1.94615e-05 0.0462434 -0.00954078 1.79158e-05 -1.80775e-05 -0.00140763
    4 ? 0.0128513 4.58248e-05 4.48221e-05 0.000127339 -0.000289767 -0.0242064 0.000105165 -0.000221029 -0.000643845 -0.000194136 -0.00102868 -1.8996e-05 -0.00165434 -0.00107047 0.0487178 2.44348e-06 -0.0183293 2.89232e-05 0.000101645 -2.39308e-05 -0.00149314

The column for *WT* identity is always *?* because the concept of a wildtype identity is not well-defined when analyzing naturally occuring sequences. The third column is the root-mean-square differential preference at that site, and the remaining columns are :math:`\Delta\pi_{r,a}` for each amino acid :math:`a` at that site.

Differential preferences absolute sum file
++++++++++++++++++++++++++++++++++++++++++++
This file has the suffix ``_diffprefsbysite_sumabs.txt`` and is created only if using *ExpCM* with ``--diffprefsbysite.txt``. It gives **half** the absolute sum the differential preferences at each site, e.g. :math:`\frac{1}{2}\sum_a \left|\Delta\pi_{r,a}\right|`. These values are sorted from largest to smallest; the possible range of values is 0 to 1. Sites with large values have large differences in preferences between natural evolution and the deep mutational scanning.

Here is an example of the first few lines::

    # The entries are half the absolute sum of the differential preferences for each site.
    # This quantity can range from 0 to 1.
    #
    # SITE  HALF_ABSOLUTE_SUM_DIFFPREFS
    61  0.36143
    5   0.276541
    15  0.255797
    58  0.244982
    12  0.237037
    26  0.217145

Dated tree file
+++++++++++++++++
This file is created only if the ``--dateseqs`` option is being used. In that case, it is the tree after least-squares dating has been applied. The file has the suffix ``_datedtree.newick``.

Most-recent common ancestor date file
+++++++++++++++++++++++++++++++++++++++++++
This file is created only if the ``--dateseqs`` option is being used. In that case, it is the time to most-recent common ancestor based on the dated tree. The file has the suffix ``_mrca_date.txt``. Its contents are very simple, for instance::

    MRCA_date = 103.135

.. include:: weblinks.txt
   
