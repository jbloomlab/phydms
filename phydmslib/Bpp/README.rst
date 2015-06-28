========================
``Bpp``
========================

This directory contains the source code for the `Bio++`_ libraries ``bpp-core``, ``bpp-seq``, and ``bpp-phyl``. These were added using ``git-subtree`` with the following commands::

    git remote add bpp-core http://biopp.univ-montp2.fr/git/bpp-core.git
    git remote add bpp-seq http://biopp.univ-montp2.fr/git/bpp-seq.git
    git remote add bpp-phyl http://biopp.univ-montp2.fr/git/bpp-phyl.git
    git subtree add --prefix=phydmslib/Bpp/bpp-core bpp-core master --squash
    git subtree add --prefix=phydmslib/Bpp/bpp-seq bpp-seq master --squash
    git subtree add --prefix=phydmslib/Bpp/bpp-phyl bpp-phyl ed4aac58830f05adfe3e192238f21849a1aa3c59 --squash

For ``bpp-phyl``, a specific commit on the ``newlik`` branch was added as that commit seems to give better performance than the subsequent patched one.



.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Bio++`: http://biopp.univ-montp2.fr/wiki/index.php/Main_Page
