.. _ExpCM:

=======================================================
**Exp**\erimentally Informed **C**\odon **M**\odels
=======================================================

.. contents::
   :depth: 2

Overview
-------------


The substitution model
------------------------

The substitution model is written in mutation-selection form. The rate of substitution :math:`P_{r,xy}` from codon :math:`x` to codon :math:`y` at site :math:`r` is

.. math::
   :label: Prxy

   P_{r,xy} = 
   Q_{xy} \times F_{r,xy} 

where :math:`Q_{xy}` is the rate of mutation from *x* to *y* and :math:`F_{r,xy}` represents the selection on this mutation. Note that the mutation terms are identical across sites, but the selection terms are site-specific.

The mutation rate terms :math:`Q_{xy}` are given by what is essentially a `HKY85`_ model. It consists of transition-transversion ratio :math:`\kappa` and four nucleotide frequencies :math:`\phi_{A}`, :math:`\phi_{C}`, :math:`\phi_{G}`, and :math:`\phi_{T}` (these nucleotide frequencies only constitute three independent parameters since they sum to one). The frequencies are the expected nucleotide composition in the **absence** of any selection on amino acids. The mutation rate terms are:

.. math::
   :label: Qxy

   Q_{xy} = \begin{cases} 
            0 & \mbox{if $x$ and $y$ differ by more than on nucleotide,} \\
            \phi_w & \mbox{if $x$ can be converted to $y$ by a transversion of a nucleotide to $w$,} \\
            \kappa \times \phi_w & \mbox{if $x$ can be converted to $y$ by a transition of a nucleotide to $w$.} \\
            \end{cases}

The selection terms :math:`F_{r,xy}` are defined in terms of the site-specific amino acid preferences. Let :math:`\pi_{r,a}` denote the preference of site :math:`r` for amino acid :math:`a`, and let :math:`\operatorname{A}\left(x\right)` denote the amino acid encoded by codon :math:`x`. Let :math:`\beta` be a *stringency parameter*, and let :math:`\omega` be a relative rate of nonsynonymous to synonymous mutations. Then the selection terms are:

.. math::
   :label: Frxy

   F_{r,xy} = 
   \begin{cases}
   1 & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
   \omega & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\pi_{r,\mathcal{A}\left(x\right)} = \pi_{r,\mathcal{A}\left(y\right)}$} \\
   \omega \times \frac{\ln\left(\left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta}\right)}{1 - \left(\left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta}\right)} & \mbox{otherwise.}
   \end{cases}

Equation :eq:`Frxy` differs from that described in `Bloom. Mol Biol Evol, 31:2753-2769`_ only by the addition of the :math:`\omega` parameter that differentiates nonsynonymous and synonymoust mutations; the equation in that reference was in turn based on the one originally derived by `Halpern and Bruno. Mol Biol Evol, 15:910-917`_ (although note that the Halpern and Bruno equation contains a key typographical error in the denominator). The *stringency parameter* :math:`\beta` has a value of greater than one if natural selection is more stringent than the experiments used to derive the preferences. The :math:`\omega` parameter is not a ratio of nonsynonymous to synonymous mutations, but rather their relative rates after accounting for the differing preferences among sites.

The stationary state of the model defined by the mutation terms alone is trivially observed to be 

.. math::
   q_x = \phi_{x_1} \times \phi_{x_2} \times \phi_{x_3}

where :math:`x_1`, :math:`x_2`, and :math:`x_3` are the nucleotides at positions 1, 2, and 3 of codon :math:`x`, and it is also trivially observed that this mutation-term only model is reversible (since direct substitution easily verifies that :math:`q_x \times Q_{xy} = q_y \times Q_{yx}`).

In `Bloom. Mol Biol Evol, 31:2753-2769`_, it is shown that the stationary state of the model defined by the selection terms alone is

.. math::
   f_x \propto \left(\pi_{r,\operatorname{A}\left(x\right)}\right)^{\beta},

and furthermore that this selection-term only model is reversible (the derivation in `Bloom. Mol Biol Evol, 31:2753-2769`_ doesn't include the :math:`\omega` term, but carrying through the same derivation with this term yields the above result).

Finally, `Bloom. Mol Biol Evol, 31:2753-2769`_ shows that given these stationary states, the overall model defined by :eq:`Prxy` is reversible and has stationary state

.. math::
   p_x = \frac{f_x \times q_x}{\sum_y f_y \times q_y} = \frac{\left(\pi_{r,\operatorname{A}\left(x\right)}\right)^{\beta} \times q_x}{\sum_y \left(\pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta} \times q_y}.



.. include:: weblinks.txt
