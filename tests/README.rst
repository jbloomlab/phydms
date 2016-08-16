==============
tests
==============

This directory contains some tests for ``phydms``. You may want to run some of these after making substantial changes to the program.

Doctests
----------

To run the simple *doctest* tests in the modules::

    python run_doctests.py


Testing ``phydms_prepalignment``
----------------------------------
Tests script for preparing codon alignments. Note that this test requires that you have `MAFFT`_ installed. Run with::

    python run_prepalignment_test.py

Testing the ``phydms`` inferences on a small NP data set
----------------------------------------------------------
To test ``phydms`` on a small data set of human influenza NPs using the preferences reported in `Doud et al (2015)`_, do::

    python run_NPtest.py

This script runs ``phydms_comprehensive`` to compare models and infer site-specific selection. It puts its results in ``./NP_test_results/`` and compares these results to the ones in ``./expected_NP_results/``. Note that this test will take a few hours to run, depending on your computer and the number of CPUs available.

Testing the ``phydms`` inferences on a Gal4 data set
----------------------------------------------------
To test ``phydms`` on a small data set of Gal4s using the preferences reported in `Kitzman et al (2014)`_, do::

    python run_Gal4test.py

This script runs ``phydms_comprehensive`` to compare models and infer site-specific selection. It puts its results in ``./Gal4_test_results/`` and compares these results to the ones in ``./expected_Gal4_results/``. Note that this test will take a few hours to run, depending on your computer and the number of CPUs available.

This data set is similar to the one used in `Bloom (2016)`_.

Testing the ``phydms`` inferences on HIV Env
----------------------------------------------
To test ``phydms`` on a small set of clade B env genes using unpublished preferences from Hugh Haddox in the Bloom lab, use::

    python run_envtest.py

This script will compare the actual results to those in ``./expected_env_results/``. It will likely take several hours to run.

Testing the ``phydms`` inferences with diversifying pressure
----------------------------------------------
To test ``phydms`` on a small set of simulated NP sequences under diversifying pressure use::

    run_divpressuretest.py
    
This script runs ``phydms`` with the flag ``--divpressure``. It puts its results in ``./divpressure_test_results`` and compares these results to the ones in ``./expected_divpressure_results``. It will likely take a few hours to run.

For more information on how the sequences were simulated see ``README.rst`` in ``./expected_divpressure_results/simulated_sequences``.

.. _`Doud et al (2015)`: https://dx.doi.org/10.1093/molbev/msv167
.. _`Kitzman et al (2014)`: http://www.nature.com/nmeth/journal/v12/n3/full/nmeth.3223.html
.. _`Bloom (2016)`: http://dx.doi.org/10.1101/037689
.. _`MAFFT`: http://mafft.cbrc.jp/alignment/software/
