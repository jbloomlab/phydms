.. _phydms_analyzeselection_prog:

=======================================
the ``phydms_analyzeselection`` program
=======================================

.. contents::

Overview
------------
``phydms_analzyeselection`` is a program that allows you to look at the distribution of site-specific selection values inferred using ``phydms`` with an *ExpCM*. It creates `Output plots`_ that show the distribution of selection values, indicates which ones are significant at a given false-discovery rate, and allows specific sites to be highlighted. 

You can run ``phydms_analyzeselection`` **after** you have used ``phydms`` (or ``phydms_comprehensive``) to create files that summarize the site-specific selection.

If you need to change the site numbers from the sequential numbering in the ``phydms`` output, use :ref:`phydms_renumber_prog` to renumber this output before passing it to ``phydms_analyzeselection``.

The `Output plots`_ created by ``phydms_analyzeselection`` are created using `matplotlib`_. You will also need command-line `latex`_ available for rendering of some of the labels, and may get an error message during the plotting if you don't have `latex`_.

See below for information on `Command-line usage`_ for ``phydms_analyzeselection``.

Command-line usage
--------------------
.. argparse::
   :module: parsearguments
   :func: PhyDMSAnalyzeSelectionParser
   :prog: phydms_analyzeselection


   outprefix
    The program will create one or more output ``PDF`` plots (one plot for each type of selection). The prefix for these files will be ``outprefix``, and the suffixes will be those described in `Output plots`_.

   selectionfiles
    You can specify an arbitrary number of files here. 

    These files must already exist, and should conform to the format of the :ref:`phydms_prog` output files that give the :math:`\omega_r`, :math:`\beta_r`, or :math:`\Delta\pi_{r,a}` values for each site (these are the ``*_omegabysite.txt``, ``*_stringencybysite.txt``, or ``*_diffprefsbysite.txt`` files, respectively). You can mix and match multiple types of file -- ``phydms_analyzeselection`` automatically detects which of these three types of selection each file gives, and puts each one in an appropriate one of the `Output plots`_.

    You can also specify multiple files for each type of selection (these can be the same gene analyzed with multiple models, or different genes), and they will be merged into the same output plot for each type of selection.

   \-\-names
    For the plots, you will typically want to assign a name to the selection in each of the ``selectionfiles``. By default, the names on the plots are simply the names of the ``selectionfiles``. However, typically you want clearer and more succinct names, so here you should list the name to go with each one of the ``selectionfiles``.

    The number of listed names must match the number of ``selectionfiles``, and each name must be unique. For multi-word names, use quotes, as in ``--names "Gal4 ExpCM" "Gal4 YNGKP"``. Note also that you can group multi-word names using the ``--groupbyname`` option.

   \-\-selectedsites
    You should use this option if you want to indicate the values for specific points on the violin plots. The value(s) given here should be the names of files that have formats like the following (which specifies that we indicate points for 7 specific sites)::

        # This table lists sites where mutations have been shown to contribute to extended-spectrum
        # resistance in TEM beta-lactamases.
        #
        # Each entry gives the site, and then and "#" followed by associated notes.
        69 # M69I, M69L, and M69V give inhibitor resistance
        130 # S130G gives inhibitor resistance
        164 # R164C, R164H, and R164S give beta-lactam resistance
        237 # A237G and A273T give beta-lactam resistance
        238 # G238S gives beta-lactam resistance
        244 # R244C, R244H, R244S give inhibitor resistance
        276 # R276D gives inhibitor resistance

    If you just give one file for ``--selectedsites``, then the sites in this file are labeled for all selections specified by ``selectionfiles``. Otherwise, you should list one file for each entry in ``selectionfiles``, and the sites in that file are labeled for each selection. If you don't want to label sites for some selections, just put ``None`` for those.

    You will get an error if the file lists site numbers for which there isn't any information in the associated file in ``selectionfiles``.

   \-\-labelselectedsites
    This option is only meaningful if you are using ``--selectedsites``. In that case, it specifies that rather than simply putting the same point marker at each site specified by ``--selectedsites``, we use a **different** marker for each point, and then place a legend that indicates the site number to the right of the plot. There is a limit to how many unique site markers are available, so you will get an error if you use this option and specify too many sites in ``--selectedsites``.

   \-\-fdr
    For the plots showing :math:`\omega_r` (the ``*_omegabysite.txt`` files) and :math:`\beta_r` (the ``*_stringencybysite.txt`` files), a blue dashed line is drawn to indicate the cutoff for sites that have value significantly different from one at this false discovery rate (computed using the `Benjamini-Hochberg`_ procedure). Use this option to specify that false discovery rate. 

    Note that the false discovery rate is computed separately for the hypothesis that the parameter is significantly > 1 and the hypothesis that the parameter is significantly < 1. 
    
    Note also that if there are no sites that are significant, the light blue line is drawn at the P-value that a site would need to have in order to be called as significant.

   \-\-maxlog10p
    For the plots showing the :math:`\omega_r` and the :math:`\beta_r` values, we show the *log10* P-values. Any *log10* P-values that have magnitude > than this are plotted as this value.

   \-\-groupbyname
    You may want to group the violin plots and show a single legend for selection values of the same type for the same gene computed with different models. If you use this option, the legend and plots are grouped for any selection files that have names (specified by ``--names``) that have the same first word. For instance, if you use ``--names "Gal4 ExpCM" "Gal4 YNGKP"`` and also specify ``--groupbyname``, then on the plot these two will be grouped together as part of the *Gal4* group. If you use this option, selection files in the same group must have the same set of selected sites specified by ``--selectedsites``.

   \-\-diffprefsline
    For the :math:`\omega_r` and :math:`\beta_r` values, we are able to determine significance using the P-values and the false discovery rate specified by ``--fdr``. But there is currently no significance testing method to determine if the absolute sum (:math:`\frac{1}{2}\sum_a \left|\Delta\pi_{r,a}\right|`) of the differential preferences plotted by this program is significantly larger than we expect.

    If you have some other method to establish when you think that the absolute sume of these differential preferences is notably large, you can use this option to draw a dotted blue line at this threshold. There are two ways of doing this:

        1) If you include some sample that you want to make your *control* sample, you can simply specify ``--diffprefsline lowestpeak``, and then the dotted blue line will be drawn at the lowest peak value for any of the samples specified in ``--selectionfiles``.

        2) If you have a numeric threshold in mind, you can specify that threshold, such as ``--diffprefsline 0.5``.


Output plots
--------------

``phydms_analyzeselection`` will generate a plot for each type of selection represented in the ``selectionfiles``. These plots will have the prefix ``outprefix`` with the following suffixes:

    * ``*_omega_violinplot.pdf`` for the :math:`\omega_r` values (see `Omega violin plot`_).

    * ``*_stringency_violinplot.pdf`` for the :math:`\beta_r` values (see `Stringency violin plot`_).

    * ``*_diffprefs_violinplot.pdf`` for the :math:`\Delta\pi_{r,a}` values (see `Differential preferences violin plot`_).

For instance, running::

    phydms_analyzeselection myplot gene1_omegabysite.txt gene2_omegabysite.txt gene1_diffprefsbysite.txt gene2_diffprefsbysite.txt

will create the output files ``myplot_omega_violinplot.pdf`` and ``myplot_diffprefsbysite.txt``. On the other hand, running::

    phydms_analyzeselection myplot gene1_omegabysite.txt gene2_omegabysite.txt 

will create just the output file ``myplot_omegaviolinplot.pdf``.

In addition to these created plots, the script prints some information to standard output about the number of sites that are significant at the FDR rates.

Omega violin plot
++++++++++++++++++++
The following command makes a simple plot ``example1_omega_violinplot.pdf`` showing the :math:`\omega_r` values::
    
    phydms_analyzeselection example1 ExpCM_omegabysite.txt YNGKP_M3_omegabysite.txt --names ExpCM YNGKP

In this plot, the y-axis shows the negative *log10* P-value that :math:`\omega_r > 1` for points with :math:`\omega_r > 1`, and the *log10* P-value that :math:`\omega_r < 1` for points with :math:`\omega_r < 1`. So points at the top of the violin plot are those for which we are more confident of diversifying selection for amino-acid change, and conversely those at the bottom are the ones for which are more confident of selection against amino-acid change. The dotted blue lines indicate the portions of the distribution for which we declare :math:`\omega_r > 1` (top line) or :math:`\omega_r < 1` (bottom line) at a false discovery rate of 0.05. When there are no significant sites, the dotted blue lines indicate the P-value that would be necessary to declare a single site significant.

.. figure:: example1_omega_violinplot.pdf
   :width: 85%
   :align: center
   :alt: If your browser does not display inline PDF images, click on this text to open in a new window.

   The plot ``example1_omega_violinplot.pdf``.

Here is a more complicated version of a similar plot that now shows to two models (grouped by gene) and labels specific points on one of them::

    phydms_analyzeselection example2 Gal4/ExpCM_omegabysite.txt Gal4/YNGKP_M3_omegabysite.txt lactamase/ExpCM_omegabysite.txt lactamase/YNGKP_M3_omegabysite.txt --names "Gal4 ExpCM" "Gal4 YNGKP" "lactamase ExpCM" "lactamase YNGKP" --selectedsites None None lactamase/resistance_sites.txt lactamase/resistance_sites.txt --labelselectedsites --groupbyname

Here is the created file ``example2_omega_violinplot.pdf``:

.. figure:: example2_omega_violinplot.pdf
   :width: 85%
   :align: center
   :alt: If your browser does not display inline PDF images, click on this text to open in a new window.

   The plot ``example2_omega_violinplot.pdf``.

Stringency violin plot
++++++++++++++++++++++++++
Similar plots to those described in `Omega violin plot`_ are created to show the :math:`\beta_r` values. Here, values are shown on the top of y-axis if there is evidence for the stringency ratio being less than one (since this indicates different selection than in the deep mutational scanning), and on the bottom if there is evidence for the stringency ratio being greater than one. For instance, the following command produces ``example3_stringency_violinplot.pdf``::

    phydms_analyzeselection example3 Gal4_stringencybysite.txt lactamase_stringencybysite.txt --names Gal4 lactamase

Here is the plot:

.. figure:: example3_stringency_violinplot.pdf
   :width: 85%
   :align: center
   :alt: If your browser does not display inline PDF images, click on this text to open in a new window.

   The plot ``example3_stringency_violinplot.pdf``.

Differential preferences violin plot
++++++++++++++++++++++++++++++++++++++
The plots showing the :math:`\beta_r` values is somewhat different. Here the y-axis shows :math:`\frac{1}{2} \sum_a \left|\Delta\pi_{r,a}\right|` for each site :math:`r`. This number can range from 0 (sites that have differential preferences of zero) to 1 (sites with maximally different differential preferences). There is no significance testing, so in the plot below we take the *Gal4* sample as our baseline, and use ``lowestpeak`` option to ``--diffprefsline`` to draw a dotted blue line at that cutoff.

Here is the command::

    phydms_analyzeselection example4 Gal4/ExpCM_diffprefsbysite.txt lactamase/ExpCM_diffprefsbysite.txt --names Gal4 lactamase --selectedsites None lactamase/resistance_sites.txt --labelselectedsites --diffprefsline lowestpeak

Here is the created plot ``example4_diffprefs_violinplot.pdf``:

.. figure:: example4_diffprefs_violinplot.pdf
   :width: 85%
   :align: center
   :alt: If your browser does not display inline PDF images, click on this text to open in a new window.

   The plot ``example4_diffprefs_violinplot.pdf``.



.. include:: weblinks.txt
