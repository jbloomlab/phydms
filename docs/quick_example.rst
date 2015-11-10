.. _quick_example:

================================
Quick example
================================

.. contents::

Input data
-----------------
To perform an analysis with `phydms`_, you need two pieces of input data:

    1) An alignment of homologous protein-coding genes.

    2) A file giving the site-specific amino-acid preferences of these genes, as might be obtained from a deep mutational scanning experiment. If you aren't sure what these preferences are, then read :ref:`ExpCM` and the references therein.

For this quick example, we will use the alignment of 37 human influenza NP sequences described in `Doud et al, Mol Biol Evol, 32:2944-2960`_. We assume that we have a FASTA file called ``humanNPs.fasta`` that contains these sequences. This file is equivalent to the human influenza sequences in Supplementary file 4 of `Doud et al, Mol Biol Evol, 32:2944-2960`_. The first few lines of this file look like this::

    >2000.01_1_A/New_York/433/2000_HOST_Human_H3N2
    ATGGCGTCCCAAGGCACCAAACGGTCTTATGAACAGATGGAAACTGATGGGGATCGCCAGAATGCAACTGAGATTAGGGCATCCGTCGGGAAGATGATTGATGGAATTGGGAGATTCTACATCCAAATGTGCACTGAACTTAAACTCAGTGATTATGAAGGGCGGTTGATCCAGAACAGCTTGACAATAGAGAAAATGGTGCTCTCTGCTTTTGATGAGAGAAGGAATAGATATCTGGAAGAACACCCCAGCGCGGGGAAAGATCCTAAGAAAACTGGAGGGCCCATATACAGGAGAGTAGATGGAAAATGGATGAGGGAACTCGTCCTTTATGACAAAGAAGAAATAAGGCGAATCTGGCGCCAAGCCAACAATGGTGAGGATGCAACAGCTGGTCTAACTCACATGATGATCTGGCATTCCAATTTGAATGATACAACATACCAGAGGACAAGAGCTCTTGTTCGAACCGGAATGGATCCCAGAATGTGCTCTCTGATGCAGGGCTCGACTCTCCCTAGAAGGTCCGGAGCTGCAGGTGCTGCAGTCAAAGGAATCGGGACAATGGTGATGGAGCTGATCAGAATGGTCAAACGGGGGATCAACGATCGAAATTTCTGGAGAGGTGAGAATGGGAGGAAAACAAGAAGTGCTTATGAGAGAATGTGCAACATTCTTAAAGGAAAATTTCAAACAGCTGCACAAAGAGCAATGGTGGATCAAGTGAGAGAAAGTCGGAACCCAGGAAATGCTGAGATCGAAGATCTCATATTTTTGGCAAGATCTGCATTGATATTGAGAGGGTCAGTTGCTCACAAATCTTGCCTACCTGCCTGTGTGTATGGACCTGCAGTATCCAGTGGGTACGACTTCGAAAAAGAGGGATATTCCTTGGTGGGAATAGACCCTTTCAAACTACTTCAAAATAGCCAAGTATACAGCCTAATCAGACCTAACGAGAATCCAGCACACAAGAGTCAGCTGGTGTGGATGGCATGCCATTCTGCTGCATTTGAAGATTTAAGATTGTTAAGCTTCATCAGAGGGACCAAAGTATCTCCGCGGGGGAAACTTTCAACTAGAGGAGTACAAATTGCTTCAAATGAGAACATGGATAATATGGGATCGAGTACTCTTGAACTGAGAAGCGGATACTGGGCCATAAGGACCAGGAGTGGAGGAAACACTAATCAACAGAGGGCCTCCGCAGGCCAAATCAGTGTGCAACCTACGTTTTCTGTACAAAGAAACCTCCCATTTGAAAAGTCAACCGTCATGGCAGCATTCACTGGAAATACGGAAGGAAGAACCTCAGACATGAGGGCAGAAATCATAAGAATGATGGAAGGTGCAAAACCAGAAGAAGTGTCGTTCCGGGGGAGGGGAGTTTTCGAGCTCTCAGACGAGAAGGCAACGAACCCGATCGTGCCCTCTTTTGACATGAGTAATGAAGGATCTTATTTCTTCGGAGACAATGCAGAAGAGTACGACAAC
    >2010.43_2_A/Managua/3244.02/2010_HOST_Human_H3N2
    ATGGCGTCCCAAGGCACCAAACGGTCTTATGAACAGATGGAAACTGATGGGGATCGCCAGAATGCAACTGAGATTAGGGCATCCGTCGGGAAGATGATTGATGGAATTGGGAGATTCTACATCCAAATGTGCACTGAACTTAAACTCAGTGATCATGAAGGGCGGTTGATCCAGAACAGCTTGACAATAGAGAAAATGGTACTCTCTGCTTTTGATGAAAGAAGGAATAAATACCTGGAAGAACACCCCAGCGCGGGGAAAGATCCCAAGAAAACTGGGGGGCCCATATACAGGAGAGTCGATGGGAAATGGATGAGGGAACTCGTCCTTTATGACAAAGAAGAAATAAGGCGAATCTGGCGCCAAGCCAACAATGGTGAGGATGCTACATCTGGTCTAACTCACATAATGATTTGGCATTCCAATTTGAATGATGCAACATACCAGAGGACAAGAGCTCTTGTTCGAACTGGAATGGATCCCAGAATGTGCTCTCTGATGCAGGGCTCGACTCTCCCTAGAAGGTCCGGAGCTGCAGGTGCTGCAGTCAAAGGAATCGGGACAATGGTGATGGAACTGATCAGAATGGTCAAACGGGGGATCAACGATCGAAATTTTTGGAGAGGTGAGAATGGGCGGAAAACAAGAAGTGCTTATGAGAGAATGTGCAACATTCTTAAAGGAAAATTTCAAACAGCTGCACAAAGAGCAATGGTGGATCAAGTGAGAGAAAGTCGGAACCCAGGAAACGCTGAGATCGAAGATCTCATATTTTTAGCAAGATCTGCATTGATATTGAGAGGATCAGTTGCTCACAAATCTTGCCTACCTGCCTGTGCGTATGGACCTGCAGTATCCAGTGGGTACGACTTCGAAAAAGAGGGATATTCCTTGGTGGGAATAGACCCTTTCAAACTACTTCAAAATAGCCAAATATACAGCTTAATCAGACCTAACGAGAATCCAGCACACAAGAGTCAGCTGGTGTGGATGGCATGCCATTCTGCTGCATTTGAAGATTTAAGATTGTTAAGCTTCATCAGAGGGACAAAAGTATCTCCTCGGGGGAAACTGTCAACTAGAGGAGTACAAATTGCTTCAAATGAGAACATGGATAATATGGGATCGAGCACTCTTGAACTGAGAAGCGGGTACTGGGCCATAAGGACCAGGAGTGGAGGAAACACTAATCAACAGAGGGCCTCCGCAGGCCAAACCAGTGTGCAACCTACGTTTTCTGTACAAAGAAACCTCCCATTTGAAAAGTCAACCATCATGGCAGCATTCACTGGAAATACGGAGGGAAGAACTTCAGACATGAGAGCAGAAATCATAAGAATGATGGAAAGTGCAAAACCAGAAGAAGTGTCATTCCGGGGGAGGGGAGTGTTCGAGCTCTCAGACGAGAAGGCAACGAACCCGATCGTGCCCTCTTTTGATATGAGTAATGAAGGATCTTATTTCTTCGGAGACAATGCAGAAGAGTACGACAAT

We will use the site-specific amino-acid preferences for NP given in `Doud et al, Mol Biol Evol, 32:2944-2960`_ (these preferences are the average of those measured in the two sequence backgrounds A/Aichi/2/1968 and A/Puerto Rico/8/1934). These preferences were obtained from deep mutational scanning data processed by `dms_tools`_ into the `preferences file format`_, and are assumed to be in a file called ``prefs.txt``. This file is equivalent to Supplementary file 3 of `Doud et al, Mol Biol Evol, 32:2944-2960`_. Here are the first few lines of this file::

    # POSITION WT SITE_ENTROPY PI_A PI_C PI_D PI_E PI_F PI_G PI_H PI_I PI_K PI_L PI_M PI_N PI_P PI_Q PI_R PI_S PI_T PI_V PI_W PI_Y
    1 M 4.16557 0.0448959 0.0461868 0.0487332 0.0461879 0.0467506 0.0411084 0.0529134 0.0270247 0.046272 0.0439419 0.172267 0.0453542 0.0454932 0.0452166 0.0447873 0.0414766 0.0379993 0.0282543 0.0464337 0.0487029
    2 A 1.99168 0.658251 0.00668612 0.0052911 0.0296401 0.00564699 0.019848 0.00767229 0.00414258 0.00780895 0.00182799 0.0367587 0.00517087 0.00934743 0.00390342 0.00596331 0.0144603 0.0138282 0.148105 0.00982134 0.00582609
    3 S 3.7758 0.0780122 0.0165141 0.00485352 0.00418567 0.155017 0.0310392 0.104461 0.0239216 0.0102993 0.107051 0.0248013 0.0594595 0.0169042 0.0205734 0.0206575 0.124961 0.101368 0.0123348 0.0251918 0.0583932

Running `phydms`_
--------------------
First, we want to see if a :ref:`ExpCM` using the preferences in ``prefs.txt`` actually describe the evolution of human influenza NP better than conventional non-site-specific substitution models. To do this, we run ``phydms_comprehensive`` using::

    phydms_comprehensive NP_analysis/ humanNPs.fasta prefs.txt --ncpus -1

Running this command creates the subdirectory ``./NP_analysis/`` which will contain the following files::

    ├── NP_analysis
    │   ├── averaged_ExpCM_prefs_diffprefsbysite.txt
    │   ├── averaged_ExpCM_prefs_diffprefsbysite_sumabs.txt
    │   ├── averaged_ExpCM_prefs.log
    │   ├── averaged_ExpCM_prefs_loglikelihood.txt
    │   ├── averaged_ExpCM_prefs_modelparams.txt
    │   ├── averaged_ExpCM_prefs_omegabysite.txt
    │   ├── averaged_ExpCM_prefs_stringencybysite.txt
    │   ├── averaged_ExpCM_prefs_tree.newick
    │   ├── ExpCM_prefs_diffprefsbysite.txt
    │   ├── ExpCM_prefs_diffprefsbysite_sumabs.txt
    │   ├── ExpCM_prefs.log
    │   ├── ExpCM_prefs_loglikelihood.txt
    │   ├── ExpCM_prefs_modelparams.txt
    │   ├── ExpCM_prefs_omegabysite.txt
    │   ├── ExpCM_prefs_stringencybysite.txt
    │   ├── ExpCM_prefs_tree.newick
    │   ├── log.log
    │   ├── modelcomparison.txt
    │   ├── YNGKP_M0.log
    │   ├── YNGKP_M0_loglikelihood.txt
    │   ├── YNGKP_M0_modelparams.txt
    │   ├── YNGKP_M0_tree.newick
    │   ├── YNGKP_M3.log
    │   ├── YNGKP_M3_loglikelihood.txt
    │   ├── YNGKP_M3_modelparams.txt
    │   ├── YNGKP_M3_omegabysite.txt
    │   └── YNGKP_M3_tree.newick

The meaning of all of these files is detailed in :ref:`phydms_prog` and :ref:`phydms_comprehensive_prog`. Briefly, files with the prefix ``ExpCM_prefs`` give the results using the *ExpCM* informed by the preferences. Files with prefix ``averaged_ExpCM_prefs`` are for a control analysis in which the preferences have been averaged across sites. The ``YNGKP_M0`` and ``YNGKP_M3`` prefixes are for analyses with the *M0* and *M3* variants of the models of `Yang, Nielsen, Goldman, and Krabbe Pederson, Genetics, 155:431-449`_.

Substitution model comparison
--------------------------------
To quickly compare the models, we can look at the ``modelcomparison.txt`` table created by ``phydms_comprehensive``. Here are the results after converting this `reStructuredText`_ table into HTML:

.. include:: example_modelcomparison.txt

Quick inspection of the `AIC`_ values in this table shows that the *ExpCM* is vastly superior to the other models. Note that the control model where the preferences are averaged across sites (*averaged_ExpCM_prefs*) is only about as good as a *YNGKP* model, showing that most of the improvement comes from capturing site-specific constraints. You can also observe that the fitted stringency parameter :math:`\beta` is greater than one for the *ExpCM*, indicating that natural evolution is more stringent than the experiment used to measure the preferences. Finally, note that :math:`\omega` (the dN/dS ratio) is quite small (:math:`\sim 0.1`) for the *YNGKP* models, but only slightly less than one for the *ExpCM*. This is because the site-specific amino-acid preferences already capture most of the contraints on the protein sequence. Overall, the results in this table are highly congruent with those in Table 1 of `Doud et al, Mol Biol Evol, 32:2944-2960`_; however, they are not entirely identical because the *ExpCM* implemented by `phydms`_ differs slightly from that used in `Doud et al, Mol Biol Evol, 32:2944-2960`_.

Identification of site-specific selection
--------------------------------------------
Sites under diversifying selection might have :math:`\omega_r` (the ratio of dN/dS) > 1. The non-site-specific models *YNGKP_M3* and *averaged_ExpCM_prefs* both identify a very large number of sites (over 100) with :math:`\omega_r < 1` at :math:`P \le 0.05`, but only one site with :math:`\omega_r > 1` (these results are in ``./NP_analysis/YNGKP_M3_omegabysite.txt`` and ``./NP_analysis/averaged_ExpCM_prefs_omegabysite.txt``). This finding is only of modest utility. Anyone with even a passing familiarity with protein evolution would be unsurprised by the fact that many sites are under some level of functional constraint (and so have a rate of nonsynonymous mutations that is less than the rate of synonymous mutations), so the fact that there are many :math:`\omega_r < 1` sites under a non-site-specific model doesn't provide much useful information. 

But when we look at the :math:`\omega_r` values obtained under the :ref:`ExpCM` model, we get much different results. This model already accounts for most of the constraints on NP's evolution, so few sites have strong evidence of :math:`\omega_r \ne 1`. Specifically, when we look at the first few lines of ``./NP_analysis/ExpCM_prefs_omegabysite.txt``, we see::

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

Now we find many fewer sites with :math:`\omega_r < 1`, but many more with :math:`\omega_r > 1`. The difference is because this comparison is now asking if sites evolve faster or slower than expected **given** some plausible model for the actual constraints operating on NP's evolution. The fast evolving sites here are therefore better candidates for being involved in some biologically important process not included in the deep mutational scanning experiments.

We can also look at the site-specific stringency values (:math:`\beta_r`) obtained under the :ref:`ExpCM` and listed in ``./NP_analysis/ExpCM_prefs_stringencybysite.txt``. Here are the first few lines::

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

We see that there are several sites for which there is very strong evidence that :math:`\beta_r < 1`. As discussed in :ref:`ExpCM`, this indicates that these sites are probably under **different** constraints in natural evolution than those measured in the lab during the deep mutational scanning. This fact is also potentially biologically interesting.

Finally, you might want to see which specific amino acids are unexpectedly favored or disfavored by natural evolution more than anticipated from the deep mutational scanning. This information is provided in the ``./NP_analysis/ExpCM_prefs_diffprefsbysite.txt`` file, which has the `dms_tools`_ `differential preferences file format`_. Here are the first few lines::

    # POSITION WT RMS_dPI dPI_A dPI_C dPI_D dPI_E dPI_F dPI_G dPI_H dPI_I dPI_K dPI_L dPI_M dPI_N dPI_P dPI_Q dPI_R dPI_S dPI_T dPI_V dPI_W dPI_Y
    1 ? 0.0054602 -0.00109088 -0.00056453 -0.00108891 -0.00125535 -0.000617395 -0.000825462 -0.00129802 -0.00107999 -0.00249272 -0.00202186 0.0236839 -0.00145916 -0.00101602 -0.00112551 -0.00222488 -0.00140373 -0.0018838 -0.000760446 -0.000350133 -0.00112508
    2 ? 0.000905382 -0.00365815 0.00010107 0.000118787 2.0868e-05 0.000105629 9.45005e-05 0.000101899 0.000182448 0.000111788 0.000119843 2.56786e-05 0.000107759 0.000111716 0.000107354 0.000119253 0.000107919 0.000269405 0.00165618 9.42552e-05 0.000101797
    3 ? 0.0117273 -0.00227315 -1.84925e-05 0.000103205 9.00981e-05 -0.0211034 -2.02e-05 -0.0042788 4.69967e-05 6.14462e-06 -0.00705743 7.46465e-06 -0.000664213 -0.00013711 -1.53571e-05 1.94615e-05 0.0462434 -0.00954078 1.79158e-05 -1.80775e-05 -0.00140763
    4 ? 0.0128513 4.58248e-05 4.48221e-05 0.000127339 -0.000289767 -0.0242064 0.000105165 -0.000221029 -0.000643845 -0.000194136 -0.00102868 -1.8996e-05 -0.00165434 -0.00107047 0.0487178 2.44348e-06 -0.0183293 2.89232e-05 0.000101645 -2.39308e-05 -0.00149314

This file is much harder to look at visually. So you might want to summarize the results using ``phydms_plotselection``. You can do this by running::

    phydms_plotselection NP_analysis/ExpCM_prefs selectionplot.pdf

This command will create the following PDF plot:

.. image:: example_phydms_plotselection_plot.pdf
   :width: 85%
   :align: center
   :alt: If your browser does not display inline PDF images, click on this text to open in a new window.

This plot shows us which sites appear to be under various forms of selection. The logo stacks, for instance, indicates that natural evolution disfavors *Asp* and favors *Gly* and *Asn* more than expected from the deep mutational scanning. We can also identify specific sites (such as 472) that are under strong diversifying selection (:math:`\omega_r > 1`) and sites (such as 442) where the stringency (:math:`\beta_r`) for the experimentally measured preferences is much lower than expected.

.. include:: weblinks.txt
