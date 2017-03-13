.. _phydms_logoplot:

=======================================
the ``phydms_logoplot`` program
=======================================

.. contents::

Overview
------------
``phydms_logoplot`` is a program that can be used to visualized amino-acid preferences or differential preferences using logo plots. 
It relies on `weblogo`_ to render the actual plots.

Running ``phydms_logoplot`` with the ``--prefs`` option allows you to visualize amino-acid preferences, which you can also re-scale by the stringency parameter (:math:`\beta`) inferred by ``phydms`` using the ``--stringency`` option.

Running ``phydms_logoplot`` with the ``--diffprefs`` option allows you to visualize the differential preferences inferred by ``phydms`` with the ``--diffprefsbysite`` option.

For both types of data, you can also overlay diversifying selection (inferred by ``phydms`` using ``--omegabysite``) via the ``--omegabysite`` option.

See below for information on `Command-line usage`_ for ``phydms_logoplot``.

Command-line usage
--------------------
.. argparse::
   :module: parsearguments
   :func: PhyDMSLogoPlotParser
   :prog: phydms_logoplot

   \-\-stringency
    Running ``phydms`` with an *ExpCM* infers a stringency parameter, :math:`\beta`. If you want your plotted preferences to be re-scaled by this stringency parameter, specify it here. A value of 1 corresponds to no re-scaling.

   \-\-nperline
    Increasing this number makes for wider and shorter plots, decreasing it makes for thinner and taller plots.

   \-\-diffprefheight
    When using ``--diffprefs``, the logo stacks will have this height in both the positive and negative direction. You will get an error if the height of a logo stack exceeds this value. Note that the maximum possible differential preference logo stack height is 1 if the preferences are completely different between the natural evolution and those provided in the site-specific amino-acid preferences file. 


Output plot
--------------
``phydms_logoplot`` prints some basic information to the screen. 
But its primary output is a PDF plot that visualizes the preferences or differential preferences.

Plotting preferences
++++++++++++++++++++++
Using the ``--prefs`` option plots amino-acid preferences.
This can either be used to visualize the preferences that you are using as input to ``phydms``, or to visualize how these preferences look after being re-scaled by the stringency parameter (:math:`\beta`) inferred by ``phydms``.

You can do this with a command like this::

    phydms_logoplot --prefs prefsfile.tsv prefs_logoplot.pdf --stringency 2.99 --nperline 72 --mapmetric charge

Here is the resulting plot:

.. image:: prefs_logoplot.pdf
   :width: 95%
   :align: center
   :alt: If your browser does not display inline PDF images, click on this text to open in a new window.

The preferences in this plot are re-scaled by the stringency parameter of :math:`\beta = 2.99`, which we got from the ``*_modelparams.txt`` file after running ``phydms``.

Note that we have set ``--mapmetric`` to color amino acids by charge. They can also be colored using other schemes.

Overlaying omega values for each site
+++++++++++++++++++++++++++++++++++++++
One use of ``phydms`` is to identify diversifying selection via the ``--omegabysite`` values. 
You may want to overlay these data onto the logoplots of either preferences or differential preferences.

Here we overlay the diversifying selection the same preferences logo plot shown above::

    phydms_logoplot --prefs prefsfile.tsv omegaoverlay_logoplot.pdf --stringency 2.99 --nperline 72 --mapmetric kd --omegabysite omegabysite.txt

where ``omegabysite.txt`` is the file produced by ``phydms``. The created plot looks like this:

.. image:: omegaoverlay_logoplot.pdf
   :width: 95%
   :align: center
   :alt: If your browser does not display inline PDF images, click on this text to open in a new window.

Note that the plot shows the *P*-value that :math:`\omega_r` is greater than or less than one (recall that the *P*-values are more meaningful than the actual :math:`\omega_r` values for this kind of analysis).
These *P*-values are **not** corrected for multiple testing (for that, look at the *Q*-values output by ``phydms``).

Note also that we have used a different ``--mapmetric``, this time coloring amino acids by hydrophobicity on the Kyte-Doolittle scale.

Plotting differential preferences
++++++++++++++++++++++++++++++++++
Using the ``--diffprefs`` option plots the differential amino-acid preferences that can be inferred by ``phydms`` using ``--diffprefsbysite``.
For each site, there is a logo stack with the height of each amino-acid single-letter code proportional to the differential preference :math:`\Delta\pi_{r,a}`. 
The y-scale is determined by ``--diffprefheight``; note that the maximum possible height in each direction is one if the amino-acid preferences are completely shifted.

Here is an example::

    phydms_logoplot --diffprefs diffprefsbysite.txt diffprefs_logplot.pdf --nperline 72

.. image:: diffprefs_logoplot.pdf
   :width: 95%
   :align: center
   :alt: If your browser does not display inline PDF images, click on this text to open in a new window.



.. include:: weblinks.txt
   
