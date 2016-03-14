.. _acknowledgments:

==================================
Acknowledgments and citations
==================================

Acknowledgments
-----------------
`phydms`_ is written by the `Bloom lab`_ (see https://github.com/jbloomlab/phydms/graphs/contributors for all contributors).

`phydms`_ is completely reliant on the `Bio++`_ libraries; thanks to `Julien Dutheil`_ and `Laurent Gueguen`_ for developing these libraries and assisting in their use.

The plotting performed by ``phydms_plotselection`` utilizes `weblogo`_.

The optional dating of nodes in the phylogenies uses the `LSD`_ source code, which was developed by the group of `Olivier Gascuel`_.

Citations
------------
There is currently no citation for `phydms`_, but the following references describe initial versions of the :ref:`ExpCM` incorporated into `phydms`_:

    * `Bloom, Mol Biol Evol, 31:1956-1978`_

    * `Bloom, Mol Biol Evol, 31:2753-2769`_

Because `phydms`_ is **completely reliant** on `Bio++`_, you should also cite these libraries. There is not yet a citation for the version of `Bio++`_ incorporated into `phydms`_, but you can cite the most recent `Bio++`_ reference:

    * `Gueguen et al, Mol Biol Evol, 30:1745-1750`_ 

If you use ``phydms_plotselection`` to create a figure for your paper, please also cite `weblogo`_ (which is used to create the logo plots):

    * `Crooks et al, Genome Research, 14:1188-1190`_

If you use the optional dating of nodes in the phylogeny, you should also cite the `LSD`_ program, which is wrapped within `phydms`_ to perform the dating. The reference for `LSD`_ is:

    * `To et al, Systematic Biology, 65:82-97`_

Finally, if you use `dms_tools`_ to process your deep mutational scanning data before feeding it into `phydms`_, please cite:

    * `Bloom, BMC Bioinformatics, 16:168`_



.. include:: weblinks.txt
