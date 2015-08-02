.. _phydms_prog:

=======================================
The ``phydms`` program
=======================================

.. contents::

Overview
------------

Command-line usage
--------------------
.. argparse::
   :module: parsearguments
   :func: PhyDMSParser
   :prog: phydms

   alignment
    This file should contain aligned DNA codon sequences. Stop codons are **not** allowed, except if a stop codon is the terminal character in all sequences, in which case they are automatically trimmed.

    The headers must be unique strings that do **not** contain any of the following: spaces, commas, colons, semicolons, parentheses, square brackets, and single or double quotes.

   tree
    This argument specifies the fixed tree topology (if **not** using ``--infertreetopology``) or an initial topology for further optimization. (if using ``--infertreetopology``). The possibilities are:

        * Specify a Newick file giving an existing tree (tip names must match sequence headers in ``alignment``).

        * Use *nj* to build a neighbor-joining tree from the nucleotide sequences using a crude identity model to compute distances. **Trees built using this option will not be very good**, so you are suggested to use this option only to get an initial tree for further optimization via ``--infertreetopology``.

        * Use *random* for a randomly chosen starting tree. If you use this option, you must also use ``--infertreetopology`` since a random tree makes no sense otherwise.

   model
    *YNGKP_M?* is a codon model from `Yang, Nielsen, Goldman, and Krabbe Pederson. Genetics, 155:431-449`_. The *?* indicates the model variant. Codon frequencies are set by *F3X4* (the product of nucleotide frequency parameters at each of the three positions; see ``--fitF3X4`` for how these 9 parameters are set). The transition-transversion ratio :math:`\kappa` is optimized by maximum likelihood. The nonsynonymous-synonymous ratio :math:`\omega` is set differently depending on the variant:
    
        - *YNGKP_M0* : a single :math:`\omega` is optimized by maximum likelihood (1 free parameter)

        - *YNGKP_M3* : three :math:`\omega` values and the weights (probabilities) of their categories optimized by maximum likelihood (5 free parameters).

        - *YNGKP_M7* : three values of :math:`\omega` between zero and one drawn from a beta distribution with the two distribution parameters optimized by maximum likelihood (2 free parameters). 

        - *YNGKP_M8* : three values of :math:`\omega` between zero and one drawn from a beta distribution with the two distribution parameters optimized by maximum likelihood, plus another value of :math:`\omega > 1` with the value and weight optimized by maximum likelihood (4 free parameters). 

    *ExpCM_<prefsfile>* is an **exp**\erimentally informed **c**\odon **m**\odel, with amino-acid preferences taken from the file ``prefsfile``. 
    The preferences file should be in the `dms_tools`_ `preferences file format`_ for **amino acids** (any stop codon preferences if present are normalized away to zero). The preferences file must specify a preference for the amino acid encoded by every site in ``alignment``, using sequential 1, 2, ... numbering.
    For information on experimentally informed codon models, see `Bloom. Mol Biol Evol, 31:2753-2769`_. 

   outprefix
    If this prefix contains a directory name, that directory is created if it does not already exist. 

    Any existing files with the names specified by ``outprefix`` are deleted at the start of the program's execution.

    The specific files created have the following suffixes appended to ``outprefix``:

        - ``.log`` : a log file recording the progress of the program

        - ``_loglikelihood.txt`` : the log likelihood after optimization

        - ``_modelparams.txt`` : the values of all optimized model parameters after optimization.

        - ``_tree.newick`` : the tree after optimization in Newick format.

   minbrlen
    Regardless of the method used to set ``tree``, all branches with lengths less than this value will be set to this value in the initial starting tree. Branches can still end up with lengths less than this after subsequent optimization of this starting tree.

   seed
    The random number seed can influence the outcome when using a non-deterministic algorithm.

   \-\-fitF3X4
    Only meaningful if ``model`` is one of the *YNGKP* models. It specifies how the 9 nucleotide frequency parameters are determined. By default, they are set to the empirically observed values in the sequence alignment. If you use ``--fitF3X4`` then these 9 parameters are fit by maximum likelihood.


.. include:: weblinks.txt
   
