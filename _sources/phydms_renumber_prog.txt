.. _phydms_renumber_prog:

=======================================
the ``phydms_renumber`` program
=======================================

.. contents::

Overview
------------
``phydms_renumber`` is a very simple program that can be used to renumber sites after analyses with ``phydms``.

This program will be useful for proteins for which the commonly used numbering scheme does **not** correspond with the sequential numbering output by ``phydms``. Specifically, ``phydms`` always outputs site numbers going 1, 2, ... for the sites in the alignment that it is provided. But for many proteins, there are common numbering schemes that might not correspond to this. If you provide this alternative numbering scheme to ``phydms_renumber`` as ``renumberfile``, it will renumber the by-site output files of ``phydms``. These are:

    * ``*_omegabysite.txt`` (if using ``--omegabysite``)

    * ``*_stringencybysite.txt`` (if using ``--stringencybysite``)

    * ``*_diffprefsbysite.txt`` and ``*_diffprefsbysite_sumabs.txt`` (if using ``--diffprefsbysite``)


Command-line usage
--------------------
.. argparse::
   :module: parsearguments
   :func: PhyDMSRenumberParser
   :prog: phydms_renumber


   renumberfile
    Site numbers do not have to be purely numerical. Here is an example file::

        1 17
        2 18
        3 18A
        4 19
        5 20

    Note that a new site number must be specified for **every** site in each file to renumber, even if that new number just matches the old one. Similarly, every site listed in this ``renumberfile`` must be present in the files to be renumbered.

    If you don't want to include a site in the output, indicate *None* for the new number, as in::

        1 None
        2 None
        3 18A
        4 19
        5 20

   outprefix
    Every file listed here should either by the full name of a file to renumber, or the prefix of such a file (we look for files with this prefix and suffixes of ``_omegabysite.txt``, ``_stringencybsite.txt``, or ``_diffprefsbysite.txt``). You are allowed to use wildcard characters (``*``) in these prefixes.

   \-\-renumberedprefix
    If ``--renumberedprefix`` is *renumbered*, then the file ``ExpCM_omegabysite.txt`` will be renumbered in the newly created file ``renumbered_ExpCM_omegabysite.txt``. The file ``mydirectory/ExpCM_omegabysite.txt`` will be renumbered into ``mydirectory/renumbered_omegabysite.txt``.



.. include:: weblinks.txt
   
