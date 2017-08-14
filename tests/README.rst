==============
tests
==============

This directory contains some tests for ``phydms``.

Each ``test_*.py`` file runs a different test, and can be executed with commands like::

    python test_expcm.py

The recommended way to run all these tests **and** the doctests in the individual modules is using `pytest`_.
To do this, first install `pytest`_ if it is not already installed.
Then navigate to the **top directory** of ``phydms`` (**not** this ``./tests/`` directory) and run `../test.bash <../test.bash>`_.


.. _`pytest`: https://docs.pytest.org/en/latest/

.. _`Doud et al (2015)`: https://dx.doi.org/10.1093/molbev/msv167
.. _`Kitzman et al (2014)`: http://www.nature.com/nmeth/journal/v12/n3/full/nmeth.3223.html
.. _`Bloom (2016)`: http://dx.doi.org/10.1101/037689
.. _`MAFFT`: http://mafft.cbrc.jp/alignment/software/
