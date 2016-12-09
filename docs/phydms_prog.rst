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
    Should contain aligned DNA codon sequences. 
    Stop codons are **not** allowed, except if a stop codon is the terminal character in all sequences, in which case they are automatically trimmed.

    No checking is done to remove identical sequences, so you may want to remove such sequences yourself if there a lot of them.

    The headers must be unique strings that do **not** contain any of the following: spaces, commas, colons, semicolons, parentheses, square brackets, and single or double quotes.

   tree
    This argument specifies the fixed tree topology used by 'phydms'. 
    Typically you would have inferred this tree using some other program such as `RAxML`_ or `codonPhyML`_.

   model
    *YNGKP_<m>* is a codon model from `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_. 
    The *<m>* indicates the model variant. 
    Codon frequencies are set empirically by *F3X4* (the product of nucleotide frequency parameters at each of the three positions).
    The transition-transversion ratio :math:`\kappa` is optimized by maximum likelihood. 
    The nonsynonymous-synonymous ratio :math:`\omega` is optimized differently depending on the variant:
     
        - *YNGKP_M0* : a single :math:`\omega` is optimized by maximum likelihood (1 free parameter)
 
        - *YNGKP_M1* : an :math:`\omega < 1` and an :math:`\omega = 1` with the weights optimized by maximum likelihood (two free parameters). This is actually the variant indicated as *M1a* in `Wong et al, Genetics, 168:1041-1051`_ rather than the original *M1* model in `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_.

        - *YNGKP_M2* : an :math:`\omega < 1`, and :math:`\omega = 1`, and an :math:`\omega > 1` with the weights optimized by maximum likelihood (4 free parameters). This is actually the variant indicated as *M2a* in `Wong et al, Genetics, 168:1041-1051`_ rather than the original *M2* model in `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_.

        - *YNGKP_M3* : three :math:`\omega` values and the weights (probabilities) of their categories optimized by maximum likelihood (5 free parameters).
 
        - *YNGKP_M7* : ``--ncats`` values of :math:`\omega` between zero and one drawn from a beta distribution with the two distribution parameters optimized by maximum likelihood (2 free parameters). 
 
        - *YNGKP_M8* : ``--ncats`` values of :math:`\omega` between zero and one drawn from a beta distribution with the two distribution parameters optimized by maximum likelihood, plus another value of :math:`\omega > 1` with the value and weight optimized by maximum likelihood (4 free parameters). 

    *ExpCM_<prefsfile>* is an :ref:`ExpCM` with amino-acid preferences taken from the file ``prefsfile``. 
    The preferences file must specify a preference for the amino acid encoded by every site in ``alignment``, using sequential 1, 2, ... numbering.
    Any stop-codon preferences specified in the file are ignored.
    The preferences file can be in several formats:

        * A comma-, tab-, or space-delimited file with the first column giving
          the site and the rest giving the preference for each amino-acid::

            site A C D E F G H I K L M N P Q R S T V W Y
            1 0.044  0.046 0.048 0.0461 0.046 0.0411 0.052 0.027 0.046 0.043 0.172 0.045 0.045 0.045 0.044 0.041 0.037 0.028 0.046 0.048
            2 0.658 0.006 0.005 0.029 0.005 0.019 0.007 0.004 0.007 0.001 0.036 0.005 0.009 0.003 0.005 0.014 0.013 0.148 0.009 0.005

        * The more complex `dms_tools`_ `preferences file format`_

   outprefix
    If this prefix contains a directory name, that directory is created if it does not already exist. 

    See `Output files`_ for a description of the files created with this prefix. Any existing files with these names are deleted at the start of the program's execution.

   \-\-brlen
    The tree topology is fixed to the that in ``tree``. 
    But the branch lengths are not necessarily fixed. 
    The options are:

      * *scale*: fix the relative length of the branches to the values in ``tree``, but scale them all by a single parameter :math:`\mu` that is optimized.

      * *optimize*: optimize all branch lengths as free parameters.

      * *fix*: fix the branch lengths exactly to their values in ``tree``. Unless you have a good reason to be sure that the lengths are already scaled correctly, use *scale* instead of *fix* if you don't want to optimize all branches.

   \-\-omegabysite
    If using a YNGKP model, then the :math:`\omega_r` value is nearly analogous that obtained using the *FEL* model described by `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_. If using and *ExpCM*, then :math:`\omega_r` has the meaning described in :ref:`ExpCM`. Essentially, we fix all other model / tree parameters and then compare a model that fits a synonymous and nonsynonymous rate to each site to a null model that only fits a synonymous rate; there is evidence for :math:`\omega_r \ne 1` if fitting both nonsynonymous and synonymous rate gives sufficiently better likelihood than fitting synonymous rate alone. See also the ``--omegabysite_fixsyn`` option.

   \-\-omegabysite_fixsyn
    This option is meaningful only if you are using ``--omegabysite``. If you use this option, then we compare a model in which we fit a nonsynonymous rate to each site to a model in which we fit nothing. The synonymous rate is not fit, and so is assumed to be equal to the overall value fit for the tree. According to `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_, in some cases this can yield greater power if there is relatively limited data. However, it comes with the risk of giving spurious results if there is substantial variation in the synonymous substitution rate among sites.  

   \-\-divpressure
    Only for *ExpCM* models. This option specifies the name of the file with the predetermined diversifying pressure for each site, :math:`\delta_{r}`. ``phydms`` will fit :math:`\omega` and :math:`\omega_{2}` as described in :ref:`ExpCM`.

    You can **not** use ``--divpressure`` and ``--fitphi`` simultaneously.

   \-\-randprefs
    Only for *ExpCM* models. 
    This option randomly reassigns the preferences among sites. 
    This can be used as a control: randomizing the preferences should make them lose their efficacy for describing evolution.

   \-\-avgprefs
    Only for *ExpCM* models. 
    This option computes an average of each preference across sites (:math:`\pi_a = \frac{1}{L} \sum_r \pi_{r,a}` where :math:`r = 1, \ldots, L`), and then uses these average preferences for all sites. 
    This can be used as a control, as it merges all the information in the preferences into a non-site-specific model.

   \-\-minbrlen
    All branches with lengths less than this value will be set to this value in the initial starting tree. 
    Branches can still end up with lengths less than this after subsequent optimization of this starting tree.

   \-\-minpref
    Adjust all amino-acid preferences in the *ExpCM* ``prefsfile`` to be at least this large.
    If this number is too close to zero, you will get runtime errors.

   \-\-seed
    The random number seed affects the outcomes when there is randomization, as when using ``--randprefs``.


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


.. include:: weblinks.txt
   
