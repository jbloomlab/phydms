Changelog
===========

1.2.0
------------
* Updated ``Bio++`` again.

1.1.0
-----------
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
