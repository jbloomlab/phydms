.. _ExpCM:

=======================================================
**Exp**\erimentally Informed **C**\odon **M**\odels
=======================================================

.. contents::
   :depth: 2

Overview
-------------
**Exp**\erimentally Informed **C**\odon **M**\odels (*ExpCM*) describe the evolution of protein-coding genes in terms of their site-specific amino-acid preferences. These models improve on conventional non-site-specific phylogenetic substitution models because they account for the different constraints at different sites of the protein encoded by the gene.

Specifically, for each gene, we assume that we know the preference :math:`\pi_{r,a}` of site :math:`r` for each amino-acid :math:`a` (we constrain :math:`1 = \sum_a \pi_{r,a}`). Typically, these preferences might be measured in deep mutational scanning experiments. For a description of how the preferences can be inferred from experimental data, see `Bloom, BMC Bioinformatics, 16:168`_. 

For examples of previous studies that have used experimentally measured site-specific amino-acid preferences in phylogenetic analyses, see:

    * `Bloom, Mol Biol Evol, 31:1956-1978`_

    * `Bloom, Mol Biol Evol, 31:2753-2769`_

    * `Thyagarajan and Bloom, eLife, 3:e03300`_

    * `Doud et al, Mol Biol Evol, 32:2944-2960`_

These experimentally informed site-specific substitution models closely parallel those used in studies that infer the site-specific information from natural sequences; see `Rodrigue and Lartillot, PNAS, 107:4629-4634`_. For a more general discussion of models of this form, see `Halpern and Bruno, Mol Biol Evol, 15:910-917`_ and `McCandlish and Stoltzfus, Quarterly Review of Biology, 89:225-252`_.

The *ExpCM* implemented by ``phydms`` are similar but not completely identical to those in aforementioned publications. The *ExpCM* in ``phydms`` has the following components:

    * A reversible model of nucleotide substitution that is assumed to be constant across sites.

    * A reversible model of amino-acid substitution that differs among sites based on the site-specific amino-acid preferences.

    * A *stringency parameter* :math:`\beta` that reflects how strongly the evolution of the gene in nature adheres to the site-specific amino-acid preferences. Under certain assumptions about the correspondence between the preferences and fitness effects of mutatons, :math:`\beta` probably has some connection to effective population size.

    * An :math:`\omega` parameter that reflects that relative rate of nonsynonymous versus synonymous substitutions **after** accounting for fact that the site-specific preferences will typically retard the overall rate of nonsynonymous substitution. This :math:`\omega` is **not** identical to the one that would be inferred using a conventional Goldman-Yang or Muse-Gaut model, since it reflects the relative rates of these two types of substitutions after accounting for the preferences. Indeed, if the preferences perfectly captured the contraints on nature evolution, we might expect :math:`\omega = 1` since there would be no further reduction in the rate of nonsynonymous substitutions after accounting for the site-specific preferences.

Below is an exact definition of the *ExpCM* model, and an explanation of how it is used to fit site-specific measures of selection.

The *ExpCM* substitution model
---------------------------------

The substitution model is in mutation-selection form. The rate of substitution :math:`P_{r,xy}` from codon :math:`x` to codon :math:`y` at site :math:`r` is

.. math::
   :label: Prxy

   P_{r,xy} = 
   Q_{xy} \times F_{r,xy} 

where :math:`Q_{xy}` is the rate of mutation from *x* to *y* and :math:`F_{r,xy}` represents the selection on this mutation. Note that the mutation terms are identical across sites, but the selection terms are site-specific.

The mutation rate terms :math:`Q_{xy}` are given by what is essentially a `HKY85`_ model. It consists of transition-transversion ratio :math:`\kappa` and four nucleotide frequencies :math:`\phi_{A}`, :math:`\phi_{C}`, :math:`\phi_{G}`, and :math:`\phi_{T}` (these nucleotide frequencies only constitute three independent parameters since they sum to one). The frequencies are the expected nucleotide composition in the **absence** of any selection on amino acids. The mutation rate terms are:

.. math::
   :label: Qxy

   Q_{xy} = \begin{cases} 
            0 & \mbox{if $x$ and $y$ differ by more than on nucleotide,} \\
            \phi_w & \mbox{if $x$ can be converted to $y$ by a transversion of a nucleotide to $w$,} \\
            \kappa \times \phi_w & \mbox{if $x$ can be converted to $y$ by a transition of a nucleotide to $w$.} \\
            \end{cases}

The selection terms :math:`F_{r,xy}` are defined in terms of the site-specific amino acid preferences. Let :math:`\pi_{r,a}` denote the preference of site :math:`r` for amino acid :math:`a`, and let :math:`\operatorname{A}\left(x\right)` denote the amino acid encoded by codon :math:`x`. Let :math:`\beta` be a *stringency parameter*, and let :math:`\omega` be a relative rate of nonsynonymous to synonymous mutations. Then the selection terms are:

.. math::
   :label: Frxy

   F_{r,xy} = 
   \begin{cases}
   1 & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
   \omega & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\pi_{r,\mathcal{A}\left(x\right)} = \pi_{r,\mathcal{A}\left(y\right)}$} \\
   \omega \times \frac{\ln\left(\left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta}\right)}{1 - \left(\left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta}\right)} & \mbox{otherwise.}
   \end{cases}

Equation :eq:`Frxy` differs from that described in `Bloom, Mol Biol Evol, 31:2753-2769`_ only by the addition of the :math:`\omega` parameter that differentiates nonsynonymous and synonymoust mutations; the equation in that reference was in turn based on the one originally derived by `Halpern and Bruno, Mol Biol Evol, 15:910-917`_ (although note that the Halpern and Bruno equation contains a key typographical error in the denominator). The *stringency parameter* :math:`\beta` has a value of greater than one if natural selection favors high-preference amino acids with greater stringency than suggested by the :math:`\pi_{r,a}` values, and has a value of less than one if natural selection favors high-preference amino acids with less stringency than suggested by the :math:`\pi_{r,a}` values.
The :math:`\omega` parameter is not a ratio of nonsynonymous to synonymous mutations, but rather their relative rates after accounting for the differing preferences among sites.

The stationary state of the model defined by the mutation terms alone is trivially observed to be 

.. math::
   q_x = \phi_{x_1} \times \phi_{x_2} \times \phi_{x_3}

where :math:`x_1`, :math:`x_2`, and :math:`x_3` are the nucleotides at positions 1, 2, and 3 of codon :math:`x`, and it is also trivially observed that this mutation-term only model is reversible (since direct substitution easily verifies that :math:`q_x \times Q_{xy} = q_y \times Q_{yx}`).

In `Bloom, Mol Biol Evol, 31:2753-2769`_, it is shown that the stationary state of the model defined by the selection terms alone is

.. math::
   f_{r,x} \propto \left(\pi_{r,\operatorname{A}\left(x\right)}\right)^{\beta},

and furthermore that this selection-term only model is reversible (the derivation in `Bloom, Mol Biol Evol, 31:2753-2769`_ doesn't include the :math:`\omega` term, but carrying through the same derivation with this term yields the above result).

Finally, `Bloom, Mol Biol Evol, 31:2753-2769`_ shows that given these stationary states, the overall model defined by :eq:`Prxy` is reversible and has stationary state

.. math::
   p_{r,x} = \frac{f_{r,x} \times q_x}{\sum_y f_y \times q_y} = \frac{\left(\pi_{r,\operatorname{A}\left(x\right)}\right)^{\beta} \times q_x}{\sum_y \left(\pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta} \times q_y}.

Therefore, assuming that the preferences :math:`\pi_{r,a}` are known *a priori* for all amino acids :math:`a` at all sites :math:`r` (e.g. the preferences have been measured in a deep mutational scanning experiment), the substitution model is completely defined by giving values to the following six parameters: :math:`\omega`, :math:`\beta`, :math:`\kappa`, :math:`\phi_A`, :math:`\phi_C`, and :math:`\phi_G`. When fitting and *ExpCM* for a gene phylogeny, ``phydms`` assumes that these six parameters are constant across all sites, and optimizes their values by maximum likelihood.


Identifying diversifying selection via site-specific :math:`\omega_r` values
------------------------------------------------------------------------------
One type of interesting selection is *diversifying selection*, where there is continual pressure for amino-acid change. Such selection might be expected to occur at sites that are targeted by adaptive immunity or subjected to some other form of selection which constantly favors changes in the protein sequence. At such sites, we expect that the relative rate of nonsynonymous substitutions will be higher than suggested by the site-specific preferences :math:`\pi_{r,a}` due to this diversifying selection.

To detect diversifying selection at specific sites within the framework of the *ExpCM* implemented in ``phydms``, we use an approach that is highly analogous the *FEL* (**f**\ixed **e**\ffects **l**\ikelihood) method described by `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_. Essentially, the tree topology, branch lengths, and all shared model parameters are fixed to their maximum-likelihood values optimized over the entire gene sequence. Then for each site :math:`r`, we fit a site-specific ratio of the rate of synonymous versus nonsynonymous substitutions while holding all holding all the other tree and model parameters constant. Effectively, this is fitting a different :math:`\omega_r` for each site, and so this analysis is indicated as ``--omegabysite`` in the ``phydms`` options.

Specifically, after fixing all of the other parameters as described above, for each site :math:`r` we re-define Equation :eq:`Frxy` as

.. math::
   :label: Frxy_omegabysite

   F_{r,xy} = 
   \begin{cases}
   \mu_r & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
   \mu_r \times \omega_r & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\pi_{r,\mathcal{A}\left(x\right)} = \pi_{r,\mathcal{A}\left(y\right)}$} \\
   \mu_r \times \omega_r \times \frac{\ln\left(\left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta}\right)}{1 - \left(\left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta}\right)} & \mbox{otherwise.}
   \end{cases}

and then fit the values for :math:`\mu_r` and :math:`\omega_r`. After this fitting, :math:`\mu_r` can be interpreted as the synonymous rate, and :math:`\mu_r \times \omega_r` as the nonsynonymous rate. The reason that we fit :math:`\mu_r` as well as :math:`\omega_r` is to model variation is synonymous rate; this can be important for the reasons described in the Discussion of `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_. If you use the ``--omegabysite_fixsyn`` option to ``phydms`` then :math:`\mu_r` is not fit, but rather is constrained to one.

The null hypothesis is that :math:`\omega_r = 1`. We compute a P-value for rejection of this null hypothesis using a :math:`\chi_1^2` test to compare the likelihood obtained when fitting both :math:`\mu_r` and :math:`\omega_r` to that obtained when fitting only :math:`\mu_r` and fixing :math:`\omega_r = 1`. See `Kosakovsky Pond and Frost, Mol Biol Evol, 22:1208-1222`_ for a justification for using a :math:`\chi_1^2` test for this type of analysis. Note that the P-values reported by ``phydms`` are **not** adjusted for multiple testing, so you will want to make such an adjustment if you are testng the hypothesis that any site has :math:`\omega_r \ne 1`. Note also that in many cases, the fitted value of :math:`\omega_r` will either be very small (e.g. close to zero) or very large (e.g. close to :math:`\infty`) -- in general, it is more informative to look for sites with small P-values and then simply look to see if :math:`\omega` is > or < 1.

Significant support for a value of :math:`\omega_r > 1` can be taken as evidence for diversifying selection beyond that expected given the constraints encapsulated in the site-specific amino-acid preferences. Significant support for a value of :math:`\omega_r < 1` can be taken as evidence for selection against amino-acid change beyond that expected given the constraints encapsulated in the site-specific amino-acid preferences. Note, however, that if the site-specific preferences don't accurately describe the real constraints, you might get :math:`\omega_r \ne 1` simply because of this fact -- so you will want to examine if sites might be subject to selection that is better described by modulating the stringency parameter :math:`\beta` or by invoking differential preferences, as described below.

Identifying differential selection via site-specific :math:`\beta_r` values
-----------------------------------------------------------------------------
Another type of interesting selection is when natural evolution prefers different amino acids than suggested by the :math:`\pi_{r,a}` values. This type of differential selection suggests discordant selection pressures between natural evolution and the process used to obtain the :math:`\pi_{r,a}` values.

One way to test for differential selection is to fit a different stringency parameter :math:`\beta_r` for each site :math:`r` after fixing the tree / branches and all other model parameters to their maximum-likelihood values obtained for the entire tree. Specifically, after fixing all of these other parameters, for each site :math:`r` we re-define Equation :eq:`Frxy` as

.. math::
   :label: Frxy_stringencybysite

   F_{r,xy} = 
   \begin{cases}
   1 & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
   \omega_r & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\pi_{r,\mathcal{A}\left(x\right)} = \pi_{r,\mathcal{A}\left(y\right)}$} \\
   \omega_r \times \frac{\ln\left(\left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta_r} / \left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta_r}\right)}{1 - \left(\left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta_r} / \left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta_r}\right)} & \mbox{otherwise.}
   \end{cases}

and then fit a value for :math:`\beta_r`. We compare the likelihood using this fitted value of :math:`\beta_r` to that obtained when using the value :math:`\beta` obtained by fitting a single :math:`\beta` over the entire gene, and compute the ratio :math:`\beta_r / \beta`. The null hypothesis is that :math:`\beta_r / \beta = 1`. We compute a P-value for rejection of this null hypothesis using a :math:`\chi_1^2` test; note that the P-values reported by ``phydms`` are **not** adjusted for multiple-hypothesis testing, so you will want to make such an adjustment if testing the hypothesis that any site has :math:`\beta_r / \beta \ne 1`. This analysis can be performed using ``phydms`` with the ``--stringencybysite`` option. Note that ``phydms`` reports the ratio :math:`\beta_r / \beta`, **not** the value of :math:`\beta_r` itself.

Significant support for :math:`\beta_r / \beta < 1` indicates that preferences don't actual describe the constraints on the gene's evolution accurately. Unfortunately, the ratio alone does not tell which particular amino acids are more or less disfavored by natural evolution, but a ratio < 1 does indicate that flatter preferences (more uniform across amino acids) do describe the evolution of the site better than the preferences encapsulated in the :math:`\pi_{r,a}` values. A value of :math:`\beta_r / \beta > 1` indicates that natural selection favors the amino acids with high :math:`\pi_{r,a}`, but that at the specific site :math:`r` this preference is more stringent than for the typical site in the gene (i.e. the preferences point in the direction of actual selection, but actual selection is more stringent at this site than suggested by the preferences).


Identifying differentially selected amino acids by fitting preferences for each site
---------------------------------------------------------------------------------------
A more complete approach is to examine each site to see the extent to which the preferences for each amino acid in nature differ from those encapsulated in the :math:`\pi_{r,a}` values. The advantage of this approach is that it can identify any form of differential selection (the approach in the previous section works best when the selection in nature is more uniform across amino acids than the :math:`\pi_{r,a}` values), and also that it can pinpoint specific amino acids that are favored or disfavored in natural evolution by an unexpected amount. The disadvantage is that ``phydms`` does not currently implement a good way to statistically test the significance of this type of differential selection, so although you can visualize and assess the selection it's hard to say that any given differential selection is significant at some specific P-value threshold.

To identify differential selection, we first fix the tree / branch lengths and all parameters of the substitution model at their maximum-likelihood values obtained across the whole tree. We then define a new set of 20 preferences for each site :math:`r` that we will denote as :math:`\hat{\pi}_{r,a}` (note that these only represent 19 parameters, as :math:`1 = \sum_a \hat{\pi}_{r,a}`. We then redefine Equation :eq:`Frxy` replacing the :math:`\left(\pi_{r,a}\right)^{\beta}` values with the :math:`\hat{\pi}_{r,a}` values as

.. math::
   :label: Frxy_diffprefs

    F_{r,xy} = 
    \begin{cases}
    1 & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
    \omega_r & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\hat{\pi}_{r,\mathcal{A}\left(x\right)} = \hat{\pi}_{r,\mathcal{A}\left(y\right)}$} \\
    \omega_r \times \frac{\ln\left(\hat{\pi}_{r,\mathcal{A}\left(y\right)} / \hat{\pi}_{r,\mathcal{A}\left(x\right)}\right)}{1 - \left(\hat{\pi}_{r,\mathcal{A}\left(x\right)} / \hat{\pi}_{r,\mathcal{A}\left(y\right)}\right)} & \mbox{otherwise.}.
    \end{cases}

Simply fitting the model defined above with these 19 :math:`\hat{\pi}_{r,a}` values will probably overfit the data since we are including 19 new parameters. We therefore regularize the parameters by defining a prior that favors :math:`\hat{\pi}_{r,a} = \left(\pi_{r,a}\right)^{\beta}`. 

Although a Dirichlet prior peaked on the preferences might seem attractive, it performs poorly in practice because the *maximum a posteriori* is very different for small preferences depending on their exact magnitude -- for instance, under a Dirichlet prior we will have very different costs of increasing the differential preference by 0.1 depending on whether the *a priori* peak estimate is :math:`10^{-3}` or :math:`10^{-4}`. This is undesirable, so instead we use a prior based on the product of inverse-quadratics.

.. math::
   :label: Pr_pi_invquad

   \Pr\left(\left\{\hat{\pi}_{r,a}\right\} \mid \left\{\pi_{r,a}\right\}, \beta\right) = \prod_a \left(\frac{1}{1 + C_1 \times \left(\hat{\pi}_{r,a} - \frac{\left(\pi_{r,a}\right)^{\beta}}{\sum_{a'} \left(\pi_{r,a'}\right)^{\beta}}\right)^2}\right)^{C_2}

or equivalently

.. math::
   :label: log_Pr_pi_invquad
    
   \log\left[\Pr\left(\left\{\hat{\pi}_{r,a}\right\} \mid \left\{\pi_{r,a}\right\}, \beta\right)\right] = - C_2 \sum_a \log\left(1 + C_1 \times \left(\hat{\pi}_{r,a} - \frac{\left(\pi_{r,a}\right)^{\beta}}{\sum_{a'} \left(\pi_{r,a'}\right)^{\beta}}\right)^2\right)

where :math:`C_1` and :math:`C_2` are concentration parameters specified via the ``--diffprefconc`` option. Larger values of these favor :math:`\hat{\pi}_{r,a}` values that more closely match the :math:`\left(\pi_{r,a}\right)^{\beta}` values (and so favor smaller values of :math:`\Delta\pi_{r,a} = \hat{\pi}_{r,a} - \frac{\left(\pi_{r,a}\right)^{\beta}}{\sum_{a'} \left(\pi_{r,a'}\right)^{\beta}}`). 

The ``phydms`` program has reasonable defaults for these concentration parameters, but you can fine tune them with ``--diffprefconc``.

Here is a plot of the log prior (Equation :eq:`log_Pr_pi_invquad`) as a function of :math:`C_1` and :math:`C_2`.

.. image:: diffprefs_prior_plot.png
   :align: center
   :width: 60%

As the plot makes clear, this prior has the desirable feature of penalizing small differential preferences proportionally more than larger ones, which is good if we think that most sites have no differential selection but some have a lot. 

To obtain the :math:`\hat{\pi}_{r,a}` values, the ``--diffprefsbysite`` option to ``phydms`` jointly optimizes the product of the prior in Equation :eq:`Pr_pi_invquad` and the likelihood (with fixed tree / branches and other model parameters) using Equation :eq:`Frxy_diffprefs`. In effect, this is the *maximum a posteriori* estimate of the :math:`\hat{\pi}_{r,a}` given the prior defined by Equation :eq:`Pr_pi_invquad`.

One might wonder why ``phydms`` takes the *maximum a posteriori* estimates rather than defining a prior with :math:`\left(\pi_{r,a}\right)^{\beta}` as the *mean* rather than the *mode*, and then sampling from the posterior to obtain the posterior mean estimates. In fact, this procedure would probably be better, but MCMC sampling of the posterior is much more computationally intensive than the *maximum a posteriori* approach.

The actual values reported by ``phydms`` are the *differential preferences*, defined as

.. math::
   :label: diffpref

   \Delta\pi_{r,a} = \hat{\pi}_{r,a} - \frac{\left(\pi_{r,a}\right)^{\beta}}{\sum_{a'} \left(\pi_{r,a'}\right)^{\beta}}.

So a differential preference of :math:`\Delta\pi_{r,a} > 0` suggests that natural evolution favors amino-acid :math:`a` at site :math:`r` more than suggested by the preferences, and a value < 0 suggests that natural evolution disfavors this amino acid. One way to summarize the total difference in preferences is the root-mean-square of the differential preferences at a site, defined as

.. math::
   :label: rms_diffpref

   RMS_{\Delta\pi_{r,a}} = \sqrt{\frac{1}{N_a} \sum_a \left(\Delta\pi_{r,a}\right)^2}.

Another way to summarize the total difference in preferences is the half the absolute sum of the differential preferences, value:

.. math::
   :label: abs_sum_diffpref

   \frac{1}{2} \sum_a \left|\Delta\pi_{r,a}\right|

The theoretical maximum of this quantity is one.

Unfortunately, ``phydms`` does not currently include any method for statistically testing (e.g. P-values) the hypothesis that any given :math:`\Delta\pi_{r,a}` is not equal to zero. So instead, you will have to keep in mind that these values are regularized (and so do not typically become substantially non-zero without some reasonable statistical evidence), and then manually inspect them for interesting trends.


.. include:: weblinks.txt
