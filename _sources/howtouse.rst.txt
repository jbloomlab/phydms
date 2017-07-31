.. _howtouse:

================================
how to use ``phydms``
================================

.. contents::

Overview
----------
The idea behind ``phydms`` is that deep mutational scanning data can inform quantitative models of protein evolution.
These models can be used to compare selection measured in the lab to that in natural evolution.

For extensive background, see the following references:

   1. `Bloom, Mol Biol Evol, 31:1956-1978`_ introduces the idea of using deep mutational scanning to build quantitative phylogenetic substitution models.
      `Bloom, Mol Biol Evol, 31:2753-2769`_ makes several extensions to this idea.

   2. `Bloom, Biology Direct, 12:1`_ describes methods to test whether sites are evolving differently in nature than expected from the deep mutational scanning.

   3. `Hilton, Doud, and Bloom, PeerJ, 2017`_ describes the ``phydms`` software package.


How to get deep mutational scanning data into a usable form 
---------------------------------------------------------------------------
The :ref:`ExpCM` used by `phydms`_ incorporate deep mutational scanning data in the form of amino-acid preferences. 
Each residue has a preference for each amino-acid, and these preferences sum to one for each residue.

Your deep mutational scanning data may or may not already be in this form. 
For instance, `dms_tools`_ directly output amino-acid preferences.
But many other methods of analyzing deep mutational scanning data give enrichment ratios or functional scores.
It is easy to convert these scores or enrichment ratios into amino-acid preferences: just normalize the enrichment ratios to sum to one at each site.
See the `Tutorial <http://htmlpreview.github.io/?https://github.com/jbloomlab/phydms/blob/master/tutorial/phydms_tutorial.html>`_ for details.

If you are missing data for some mutations, you have to interpolate them somehow.
If you are only missing a few mutations, a reasonable approach is to estimate their effect as equal to the average for all other mutations.
You can't simply leave some estimates out -- any mutation can happen in evolution, so an evolutionary model has to include an estimate for each mutation's effects.

You can **not** simply estimate amino-acid preferences as the frequencies in a natural sequence alignment and use them in `phydms`_. 
That approach is circular: the preferences need to come from data independent of the natural sequences.
A deep mutational scanning experiment is independent, a natural alignment is probably not. 
If you want to try to get the preferences out of the alignment, consider Bayesian approaches like the ones `here <https://academic.oup.com/mbe/article/34/1/204/2656188/Detecting-Adaptation-in-Protein-Coding-Genes-Using>`_ `or here <https://academic.oup.com/bioinformatics/article/30/7/1020/234902/Site-heterogeneous-mutation-selection-models>`_; these  
Bayesian approaches `don't have the same problem with overfitting <http://bayesiancook.blogspot.com/2014/01/the-myth-of-over-parameterization.html>`_, although they are vastly more computationally expensive and have other shortcomings.


Tutorial
-----------
The best way to learn to use `phydms`_ is to look at the tutorials.
`You can access the tutorials by clicking here <http://htmlpreview.github.io/?https://github.com/jbloomlab/phydms/blob/master/tutorial/phydms_tutorial.html>`_.

.. include:: weblinks.txt
