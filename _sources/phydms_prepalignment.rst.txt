.. _phydms_prepalignment_prog:

=======================================
the ``phydms_prepalignment`` program
=======================================

.. contents::

Overview
------------
``phydms_prepalignment`` is a program that is useful for preparing your sequence alignment for input to ``phydms``.

It can help you make **codon**-level alignments from protein alignments of translated coding sequences. It also performs basic checks on sequences, and can remove overly diverged or similar sequences. The basic notion is that you have some reference sequence (presumably what you used for the deep mutational scanning), and want to get homologs aligned to that for analysis.

See below for information on `Command-line usage`_ and `Examples`_.

Command-line usage
--------------------
.. argparse::
   :module: parsearguments
   :func: PhyDMSPrepAlignmentParser
   :prog: phydms_prepalignment


Examples
----------
Let's say we have some unaligned homologous protein-coding sequences in ``SUMO1_orthologs.fasta``. Let's say we want to use as our reference sequence the human ortholog, which has the header ``>ENSP00000376077_Hsap/1-303`` (note the *Hsap*, indicating *Homo sapiens*). We can build a codon alignment with::

    phydms_prepalignment SUMO1_orthologs.fasta SUMO1_alignment.fasta Hsap --minidentity 0.7 --minuniqueness 2

This creates the alignment ``SUMO1_alignment.fasta`` in which all gaps relative to the *Homo sapiens* orthologs are stripped, all sequences have at least 70% identity to the *Homo sapiens* ortholog, and all sequences have at least two protein differences from other retained homologs.

In addition to creating the alignment, the plot ``SUMO1_alignment.pdf`` is crated which shows the nucleotide and protein divergence from the reference *Homo sapiens* sequence of all homologs that passed the filters related to no premature stops, divisible by 3, and no stop codons. Here is that plot:

.. figure:: SUMO1_alignment.jpg
   :width: 40%
   :align: center
   :alt: SUMO1_alignment.jpg

   A ``jpg`` version of the plot ``SUMO1_alignment.pdf``.

Sometimes there may be certain sequences you do **not** want to purge even if they lack uniqueness or identity to the reference sequence. For instance, maybe you want to keep the chimpanzee and gorilla homologs (which have *Panu* and *Ggor* in their headers) even though they are less than two amino-acid mutations from the human ortholog. In this case, use ``--keepseqs``::

    phydms_prepalignment SUMO1_orthologs.fasta SUMO1_alignment.fasta Hsap --minidentity 0.7 --minuniqueness 2 --keepseqs Panu Ggor

If there are many such sequences, you can also list them in a file, such as ``keepseqsfile.txt``::

    Panu
    Ggor

and then specify this file to ``--keepseqs`` as::

    phydms_prepalignment SUMO1_orthologs.fasta SUMO1_alignment.fasta Hsap --minidentity 0.7 --minuniqueness 2 --keepseqs keepseqsfile.txt

If there are sequences that you specifically want to exclude, you can do something similar with ``--purgseqs``.


.. include:: weblinks.txt
   
