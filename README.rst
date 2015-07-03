========================
``phydms``
========================
``phydms`` enables **phy**\logenetic analyses using **d**\eep **m**\utational **s**\canning data to inform the substitution models.

Written by `Jesse Bloom`_.

Installation from source
------------------------
To install from source, if you are superuser do the following to install globally::

    python setup.py install

To install locally (to ``~/.local/``) use::

    python setup.py install --user

Bio++ code 
----------------------------------
The `Bio++`_ libraries are used for the phylogenetic likelihood calculations. 

To avoid requiring separate installation of the `Bio++`_ libraries, their source code is included in the source for ``phydms`` and they are compiled and statically linked to ``phydms``. Specifically, the ``./phydmslib/Bpp/`` subdirectory holds the ``bpp-core``, ``bpp-seq``, and ``bpp-phyl`` libraries, which were added to the repository using ``git subtree``. These libraries are compiled by ``cythonize`` as specified in the ``setup.py`` file. 

The need to compile and statically link the `Bio++`_ libraries increases the time to build ``phydms``. The ``setup.py`` script includes an option that should allow you to instead dynamically link to pre-compiled `Bio++`_ libraries. For this to work, you will need to have pre-installed the correct versions of ``bpp-core``, ``bpp-seq``, and ``bpp-phyl`` in the location expected by ``setup.py``. Then run::

    python setup.py install --user --dynamically-link-bpp

or::

    python setup.py install --dynamically-link-bpp

Note that this ``--dynamically-link-bpp`` option may lead to crashes during installation or execution if you don't have the right versions of `Bio++`_. Therefore, only use this option if you know what you're doing.



.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Bio++`: http://biopp.univ-montp2.fr/wiki/index.php/Main_Page
