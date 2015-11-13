.. _phydms_plotselection_prog:

=======================================
the ``phydms_plotselection`` program
=======================================

.. contents::

Overview
------------
``phydms_plotselection`` is a program that visually displays the site-specific selection inferred using ``phydms`` with an *ExpCM*. Specifically, it creates an `Output plot`_ that visualizes the following forms of selection:

    - *differential selection* : this is the differential selection between natural evolution and the site-specific preferences for each amino acid at each site, as represented by the :math:`\Delta\pi_{r,a}` values described in :ref:`ExpCM`.

    - *stringency by site* : this is the stringency with which natural evolution adheres to the site-specific preferences at each site, as represented by the :math:`\beta_r` values described in :ref:`ExpCM`.

    - *omega by site* : this is the diversifying selection for nonsynonymous over synonymous change at each site, as represented by the :math:`\omega_r` values described in :ref:`ExpCM`.

You can run ``phydms_plotselection`` **after** you have used ``phydms`` (or ``phydms_comprehensive``) to infer the site-specific selection -- ``phydms_plotselection`` requires that the output ``phydms`` files describing the site-specific selection already exist. 

If you need to change the site numbers from the sequential numbering in the ``phydms`` output, use :ref:`phydms_renumber_prog` to renumber this output before passing it to ``phydms_plotselection``.

The `Output plot`_ created by ``phydms_plotselection`` is created using `matplotlib`_ and the `weblogolib`_ that comes with `weblogo`_. 
See below for information on `Command-line usage`_ for ``phydms_plotselection``.

Command-line usage
--------------------
.. argparse::
   :module: parsearguments
   :func: PhyDMSPlotSelectionParser
   :prog: phydms_plotselection


   plotfile
    The file created with this name is described in `Output plot`_.

   \-\-nperline
    Increasing this number makes for wider and shorter plots, decreasing it makes for thinner and taller plots.

   \-\-diffprefheight
    The logo stacks in the `Output plot`_ will have this height in both the positive and negative direction. You will get an error if the height of a logo stack exceeds this value. Note that the maximum possible differential preference logo stack height is 1 if the preferences are completely different between the natural evolution and those provided in the site-specific amino-acid preferences file. If you use ``--updiffprefheight`` then the logo stack height will be automatically increased (but **never** decreased) from the value specified here if necessary.

   \-\-updiffprefheight
    Use this option if you want to automatically increase the logo stack height if the differential preferences at any site exceed ``--diffprefheight``.

   \-\-minP
    The overlays showing the stringency-by-site ratio and omega-by-site are colored by the *P*-value with which we can reject the null hypothesis that these quantities are equal to one. This coloring will span *P* values from 1 to ``--minP``; any values less than ``--minP`` will be colored with the same color as ``--minP``.


Output plot
--------------
``phydms_plotselection`` prints some basic information to the screen. But its primary output is a PDF plot that visualizes the differential selection. For each site, there is a logo stack with the height of each amino-acid single-letter code proportional to :math:`\Delta\pi_{r,a}`. In addition, if you are using ``--omegabysite`` or ``--stringencybysite`` there will be overlay bars that use colors to show the *P*-value with which we can reject the null hypothesis that :math:`\omega_r` and :math:`\beta_r / \beta` are not equal to one. These color bars indicate whether the fitted value of :math:`\omega_r` and :math:`\beta_r / \beta` are greater or less than one.

Here is an example of a created plot:

.. image:: example_phydms_plotselection_plot.pdf
   :width: 85%
   :align: center
   :alt: If your browser does not display inline PDF images, click on this text to open in a new window.

.. include:: weblinks.txt
   
