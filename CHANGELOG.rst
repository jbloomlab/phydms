Changelog
===========

1.2.3
----------
* Updated docs

* Handle negative MRCAs for ``--dateseqs`` without error

* Added option ``--colormap`` to ``phydms_plotselection``.

* Removed requirements that input sequences be unique in ``alignment``.

1.2.2
--------
* Fixed bug in handling *YNGKP_M7* model in output of ``phydms_comprehensive``

* Added option ``--ncats`` option to set the number of categories for the distributed *omega* in the *YNGKP_M7* and *YNGKP_M8* models.

* Fixed bug in ``--ngammarates`` that caused only one rate to be used.

1.2.1
----------
* Eliminate negative branch lengths to avoid bug in ``biopython`` (<= version 1.66) at parsing them.

1.2.0
------------
* Updated ``Bio++`` again.

* Added ``--dateseqs`` option ``phydms`` and ``phydms_comprehensive``.

* Fixed bug that was making *YNGKP* models always use ``--fitF3X4`` option.

* Added *YNGKP_M1* and *YNKGP_M2* models.

* Changed meaning of ``--yngkp`` option to ``phydms_comprehensive`` to allow multiple or no models beyond *YNGKP_M0*.

* Added ``--gammarates`` option to ``phydms`` and ``phydms_comprehensive``.

* Added env test to ``./tests/``

* Changed source code and docs url from *jbloom* to *jbloomlab* account on ``GitHub``.

1.1.0
-----------
* Updated to newest versions of ``Bio++``

* Added ``--yngkp`` option to ``phydms_comprehensive``, and changed default from *YNGKP_M3* to *YNGKP_M8*.

* Added ``--useLog`` as option to ``phydms`` and ``phydms_comprehensive``, and made differential preferences automatically try logarithms when it encounters problems.

* Added ``--avgrandcontrol`` option to ``phydms_comprehensive``.

* Disallowed identical sequences in ``alignment``.

* Updated to newer versions of ``Bio++``

* Added output for site likelihoods to log file under the fixed stringency parameter when using ``--stringencybysite``.

* Fixed ``--no_optimize`` method for *YNGKP_M3*, *YGNKP_M7*, and *YNGKP_M8* by making program use the old likelihood method for these models.


1.0.2
--------
* Fixed bug in file checks for ``--no_optimize`` option to ``phydms``

* Add ``--no_avgprefs`` option for ``phydms_comprehensive``

1.0.1
--------
* Included ``__*`` files in ``Bpp`` in MANIFEST for proper ``PyPI`` / ``pip`` installation

1.0.0
--------
Initial release
