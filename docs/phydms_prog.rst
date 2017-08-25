.. _phydms_prog:

=======================================
the ``phydms`` program
=======================================

.. contents::

Overview
------------
``phydms`` can perform phylogenetic analyses with :ref:`ExpCM` as well as with some standard non-site-specific substitution models (variants of the *YNGKP* models described in `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_). In addition, ``phydms`` can be used to detect site-specific diversifying or differential selection using the :ref:`ExpCM` as described in `Bloom, Biology Direct, 12:1`_.

See below for information on the `Command-line usage`_ and `Output files`_.

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
    Typically you infer this tree using a program such as `RAxML`_ or `codonPhyML`_.
    The topology is fixed to that in ``tree``. How the branch lengths are handled
    during optimization depends on ``--brlen``.
    The branch lengths in ``tree`` are assumed to be in codon substitutions
    per site.

   model
    *YNGKP_<m>* is a *non-site-specific* Goldman-Yang style codon model from `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_.
    The *<m>* indicates the model variant.
    Codon frequencies are set empirically using the *CF3X4* method described by `Pond et al, PLoS One, 5:e11230`_ (this sets the codon frequencies to the the product of nucleotide frequency parameters at each of the three positions after correcting for stop codons).
    The transition-transversion ratio :math:`\kappa` is optimized by maximum likelihood.
    The nonsynonymous-synonymous ratio :math:`\omega` is optimized differently depending on the variant:

        - *YNGKP_M0* : a single :math:`\omega` is optimized by maximum likelihood.

        - *YNGKP_M5* : :math:`\omega` drawn from a gamma distribution with the two parameters optimized by maximum likelihood. The distribution has ``--ncats`` categories.

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

    Importantly, the amino-acid preferences must be obtained **independently** from the sequence
    alignment being analyzed. A deep mutational scanning experiment is an independent means
    of obtaining the preferences, but estimating them from the amino-acid frequencies in the alignment of homologs is not a valid approach as you are then estimating the preferences from the same sequences that you are subsequently analyzing.

   outprefix
    If this prefix contains a directory name, that directory is created if it does not already exist.

    See `Output files`_ for a description of the files created with this prefix. Any existing files with these names are deleted at the start of the program's execution.

   \-\-brlen
    The tree topology is fixed to that in ``tree``.
    But there are several options for how the branch lengths are handled:

      * *scale*: fix the relative length of the branches to the values in ``tree``, but scale them all by a single parameter :math:`\mu` that is optimized. This approach will be faster.

      * *optimize*: optimize all branch lengths as free parameters. This approach will be more accurate, but slower.

   \-\-gammaomega
    Use this option for a `model` of *ExpCM* if you want :math:`\omega` to be drawn from ``--ncats`` gamma-distributed categories, similar to a *YNGKP_M5* model.
    This option adds one extra free parameter, as we are now optimizing two parameters for the :math:`\omega` distribution. It will also increase run-time by about 5-fold if you use the default for ``--ncats``.

    You do **not** use this option with *YNGKP* models; for those, simply set ``model`` to *YNGKP_M5* to get gamma-distributed :math:`\omega`.

   \-\-ncats
    Determines the number of categories when using ``--gammaomega`` and `model` of *YNGKP_M5*. More categories leads to longer run-time, values of 4-5 are usually adequate.

   \-\-fitphi
    This option is not typically recommended. It will typically lead to only very slight improvements in log likelihood at substantial computational cost.

   \-\-omegabysite
    If using a YNGKP model, then the :math:`\omega_r` value is nearly analogous that obtained using the *FEL* model described by `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_. If using and *ExpCM*, then :math:`\omega_r` has the meaning described in :ref:`ExpCM`. Essentially, we fix all other model / tree parameters and then compare a model that fits a synonymous and nonsynonymous rate to each site to a null model that only fits a synonymous rate; there is evidence for :math:`\omega_r \ne 1` if fitting both nonsynonymous and synonymous rate gives sufficiently better likelihood than fitting synonymous rate alone. See also the ``--omegabysite_fixsyn`` option.

   \-\-omegabysite_fixsyn
    This option is meaningful only if you are using ``--omegabysite``. If you use this option, then we compare a model in which we fit a nonsynonymous rate to each site to a model in which we fit nothing. The synonymous rate is not fit, and so is assumed to be equal to the overall value fit for the tree. According to `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_, in some cases this can yield greater power if there is relatively limited data. However, it comes with the risk of giving spurious results if there is substantial variation in the synonymous substitution rate among sites.

   \-\-diffprefsbysite
    This option can only be used with *ExpCM* models, **not** with *YNGKP* models.

    The fitting of site-specific differential selection is designed to identify sites that prefer different amino acids in the natural alignment than expected from the deep mutational scanning data.
    This approach is described in `Bloom, Biology Direct, 12:1`_.

    The prior that is used is specified by ``--diffprefsprior``.

   \-\-diffprefsprior
    How do we regularize the differential preferences when using ``--diffprefsbysite``?

    To use the inverse-quadratic prior in `Bloom, Biology Direct, 12:1`_, specify
    the string `invquadratic,C1,C2` where `C1` and `C2` are numbers > 0 that specify
    how concentrated the prior is.

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

   \-\-fitprefsmethod
    As described in :ref:`implementation`, there are two different internal parameter representations that can be used when optimizing the preferences as free parameters.
    Both in principle should give the same result, but in practice one method might be more efficient than the other or work better in rare cases where the optimization runs into problems.

   \-\-seed
    The random number seed affects the outcomes when there is randomization, as when using ``--randprefs``.

   \-\-ncpus
    Using multiples CPUs accelerates the by-site fitting (e.g., ``--omegabysite``) as multiple sites can be analyzed at once.
    However, it does not accelerate the initial optimization of the model and tree, as this calculation is not currently parallelized.


Output files
--------------
Running ``phydms`` will create the following output files, all with the prefix specified by ``outprefix``.

Log file
+++++++++++++
This is a text file with the suffix ``_log.log`` that records information about the program's progress.

Log likelihood file
++++++++++++++++++++++
This file has the suffix ``_loglikelihood.txt``. It simply gives the optimized log likelihood. Here is an example of a file's contents::

    log likelihood = -4415.24

Model parameters
++++++++++++++++++++
This file has the suffix ``_modelparams.txt``.

Here is an example of the contents for an *ExpCM* model.
It gives the values of the parameters described in :ref:`ExpCM`; note that :math:`\mu` is **not** included as it is confounded with the branch lengths::

    beta = 2.99307
    kappa = 6.31364
    omega = 0.780782
    phiA = 0.372278
    phiC = 0.196416
    phiG = 0.224367

Here is an example of the contents for a *YNGKP_M0* model.
In this file, `phi0A` is the corrected empirical frequency of `A` at the first codon site, `phi0A` is the frequency of `A` at the second codon site, etc.::

    kappa = 5.78484
    omega = 0.11308
    phi0A = 0.340656
    phi0C = 0.153395
    phi0G = 0.308294
    phi1A = 0.320503
    phi1C = 0.206698
    phi1G = 0.236844
    phi2A = 0.331485
    phi2C = 0.202991
    phi2G = 0.243997

If you use a model with a gamma-distributed :math:`\omega` (i.e., the ``--gammarates`` option for an *ExpCM*, or the *YNGKP_M5* model) or :math:`\beta`, rather than have a single value for the parameter, there are instead two parameters that determine the gamma distribution.
For a gamma-distributed :math:`\omega`, these are the shape parameter :math:`\alpha_{\omega}` (denoted *alpha_omega*) and the inverse scale parameter :math:`\beta_{\omega}` (denoted by *beta_omega*).
The mean and variance of the omega distribution are :math:`\alpha_{\omega}/ \beta_{\omega}` and :math:`\alpha_{\omega} / \left(\beta_{\omega}\right)^2`, respectively.
To get the exact values, use the :ref:`api` to call ``phydmslib.models.DiscreteGamma(alpha_omega, beta_omega, ncats)`` where *ncats* is the value set by ``--ncats``. Here is an example of the model parameter file contents for an *ExpCM* with ``--gammaomega``::

    alpha_omega = 0.835183
    beta = 3.01549
    beta_omega = 0.976053
    kappa = 6.3424
    phiA = 0.372475
    phiC = 0.196471
    phiG = 0.22418

The shape and scale parameters for gamma-distributed :math:`\beta` are :math:`\alpha_{\beta}` and :math:`\beta_{\beta}`, respectively. 

Tree file
++++++++++++
This file has the suffix ``_tree.newick``, and gives the optimized tree in Newick format.
The branch lengths are in units of codon substitutions per site.

Site-specific omega file
+++++++++++++++++++++++++++
This file has the suffix ``_omegabysite.txt``, and is created only if using the ``--omegabysite`` option. This file gives the :math:`\omega_r` values optimized for each site. In the case of a *YNGKP* model, these are site-specific dN/dS ratios that should be essentially analogous to those obtained under the *FEL* model described by `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_. In the case of an *ExpCM* model, these values indicate diversifying selection for nonsynonymous amino-acid change as described in `Bloom, Biology Direct, 12:1`_ (see also :ref:`ExpCM`).

Here is an example of the first few lines of a file. The entries are tab separated::

    site    omega   P   dLnL    Q
    213 0.000   0.000216    6.841   0.108
    55  0.000   0.000649    5.815   0.162
    354 0.000   0.00166 4.948   0.177
    298 0.000   0.00181 4.865   0.177

The *P*-values are for rejection of the null hypothesis that :math:`\omega_r = 1` (calculated using a :math:`\chi_1^2` test).
These *P*-values are **not** corrected for multiple hypothesis testing; so you should instead look at the *Q*-values, which give the false discovery rate if this site and all above it are considered to have :math:`\omega \ne 1`.
The sites are sorted in the file by *P*-value.

For this type of test, the quantity of interest is typically not so much the exact value of :math:`\omega_r`, but rather the *P*-value (or *Q*-value) that :math:`\omega_r \ne 1`.
That *P*-value (or *Q*-value) is what indicates the credence that you should lend to the idea that a site is under diversifying selection if :math:`\omega_r > 1`.

Site-specific differential selection
+++++++++++++++++++++++++++++++++++++++
This file has the suffix ``_diffprefsbysite.txt``, and is created only if using the ``--diffprefsbysite`` option to an *ExpCM*. These are the differential preferences (the :math:`\Delta\pi_{r,a}` values) described in `Bloom, Biology Direct, 12:1`_ and :ref:`ExpCM`.

Here is an example of the first few lines of a file. The entries are tab separated::

    site    dpi_A   dpi_C   dpi_D   dpi_E   dpi_F   dpi_G   dpi_H   dpi_I   dpi_K   dpi_L   dpi_M   dpi_N   dpi_P   dpi_Q   dpi_R   dpi_S   dpi_T   dpi_V   dpi_W   dpi_Y   half_sum_abs_dpi
    375 -0.0004 -0.0040 0.0995  -0.8408 -0.0014 0.4839  -0.0029 -0.0003 -0.0062 -0.0333 -0.0092 -0.0039 -0.0003 -0.0340 -0.0000 -0.0036 -0.0073 0.3659  -0.0013 -0.0003 0.9492
    421 -0.0107 -0.0062 -0.0092 0.1521  -0.0009 -0.0004 -0.0078 -0.0124 -0.0003 -0.0168 -0.0036 -0.0271 -0.0000 -0.0019 -0.0060 -0.0259 -0.0022 -0.0116 -0.0003 -0.0085 0.1521
    283 -0.0126 -0.0124 -0.0001 -0.0031 -0.0052 -0.0129 -0.0131 -0.0047 -0.0014 -0.0085 -0.0077 -0.0012 0.1502  -0.0123 -0.0004 -0.0237 -0.0164 -0.0089 -0.0054 -0.0000 0.1502
    294 -0.0078 -0.0087 -0.0112 0.1458  -0.0084 -0.0079 -0.0098 -0.0001 -0.0064 -0.0028 -0.0079 -0.0134 -0.0000 -0.0126 -0.0104 -0.0156 -0.0110 -0.0115 -0.0001 -0.0002 0.1458
    127 -0.0088 -0.0006 -0.0010 0.1423  -0.0021 -0.0179 -0.0059 -0.0096 -0.0208 -0.0100 -0.0021 -0.0095 -0.0007 -0.0066 -0.0073 -0.0114 -0.0146 -0.0075 -0.0010 -0.0049 0.1423
    289 -0.0079 -0.0127 -0.0005 -0.0002 -0.0228 -0.0005 -0.0154 -0.0156 -0.0033 -0.0167 -0.0113 -0.0034 -0.0004 -0.0004 -0.0094 -0.0020 -0.0028 -0.0133 -0.0006 0.1391  0.1391

The first column gives the site numbers, subsequent columns give the differential preference (:math:`\Delta\pi_{r,a}`) for each amino acid.
The last column gives the half absolute sum of the differential preferences, :math:`\sum_a |\Delta\pi_{r,a}|`, at each site. This quantity can range from zero to one.
The sites are sorted with the highest half absolute sum differential preference first.


.. include:: weblinks.txt
