========================
``phydms``
========================
``phydms`` enables **phy**\logenetic analyses using **d**\eep **m**\utational **s**\canning data to inform the substitution models.

Written by `Jesse Bloom`_.

`Bio++`_ code 
----------------------------------
The `Bio++`_ libraries are used for the phylogenetic likelihood calculations. In order to avoid requiring separate installation of these libraries, their source code is included in the source code for ``phydms``. Specifically, the ``./phydmslib/Bpp/`` subdirectory holds the ``bpp-core``, ``bpp-seq``, and ``bpp-phyl`` libraries, which were added to the repository using ``git subtree``. These libraries are compiled by ``cythonize`` as specified in the ``setup.py`` file. This compilation substantially increases the time to build and install ``phydms``. If you want to link to pre-compiled `Bio++`_ libraries instead, you will need to modify the definition of the extensions for ``cythonize`` in ``setup.py``.


.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Bio++`: http://biopp.univ-montp2.fr/wiki/index.php/Main_Page
