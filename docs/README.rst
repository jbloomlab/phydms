===========================
Documentation
===========================

This subdirectory contains the `reStructuredText`_ documentation that can be built with `sphinx`_.

Building the documentation
-----------------------------

To build the documentation, you will need to have installed:

    * `sphinx`_ 
    
    * `sphinx-argparse`_ 

    * `sphinx-apidoc`_ 

Then simply type::

    make html

and the HTML documentation will be installed in ``./_build/html/``.

Notes
--------

Note that the configuration automatically created by ``sphinx-quickstart`` has been modified in the following ways:

    * ``conf.py`` has been modified to read the version and package information from ``../phydmslib/_metadata.py``

    * ``Makefile`` has been modified to automatically run `sphinx-apidoc`_ when invoked with ``make html``.


.. _`reStructuredText`: http://docutils.sourceforge.net/docs/user/rst/quickref.html
.. _`sphinx`: http://sphinx-doc.org/
.. _`sphinx-apidoc`: http://sphinx-doc.org/man/sphinx-apidoc.html
.. _`sphinx-argparse`: http://sphinx-argparse.readthedocs.org
