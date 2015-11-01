==================
LSDExtensions
==================

This directory contains ``c++`` files that enable `LSD`_ to be extended so that it can be used by ``phydms`` via a `cython`_ wrapper.

``lsd.cpp`` is the same file copied from the `LSD`_ source with a few simple modification that the program name is changed from ``main`` to ``lsd``. Essentially, this is the executable for the `LSD`_ program with slight re-purposing to enable better `cython`_ wrapping.

``lsd.h`` is a header for ``lsd.cpp``.

.. _`LSD`: http://www.atgc-montpellier.fr/LSD/
.. _`cython`: http://cython.org/
