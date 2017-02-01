.. _YNGKP_M0:

=======================================================
YNGKP_M0
=======================================================
.. contents::
   :depth: 2

Substitution model parameters and derivatives
---------------------------------------------------------
We begin by considering the basic *YNGKP_M0* substitution model defined in Goldman and Yang, 1994.
:math:`P_{xy}` is the substitution rate from codon `x` to codon `y` and is defined by

.. math::
  :label: Pxy_M0

  P_{xy} =
  \begin{cases}
  0 & \mbox{if $x$ and $y$ differ by more than one nucleotide,}\\
  \mu \omega \Phi_{y} & \mbox{if $x$ is converted to $y$ by a single-nucleotide transversion,} \\
  \kappa \mu \omega \Phi_{y} & \mbox{if $x$ is converted to $y$ by a single-nucleotide transition,} \\
  -\sum\limits_{z \ne x} P_{xz} & \mbox{if $x = y$.}
  \end{cases}

where :math:`\kappa` is the transition-transversion ratio, :math:`\Phi_y` is the empirical frequency of
codon :math:`y`, :math:`\omega` is the gene-wide rate of non-synonymous change, and :math:`\mu` is the substitution rate.
The free parameters are :math:`\kappa` and :math:`\omega`.

The derivatives are:

.. math::
   :label: dPxy_dkappa_M0

   \frac{\partial P_{xy}}{\partial \kappa}
   &=&
   \begin{cases}
   \frac{P_{xy}}{\kappa} & \mbox{if $x$ is converted to $y$ by a transition of a nucleotide to $w$,} \\
   0 & \mbox{if $x$ and $y$ differ by something other than a single transition,} \\
   -\sum\limits_{z \ne x} \frac{\partial P_{xz}}{\partial \kappa} & \mbox{if $x = y$.}
   \end{cases}

.. math::
      :label: dPxy_domega_M0

      \frac{\partial P_{xy}}{\partial \omega} =
      \begin{cases}
      0 & \mbox{if $\operatorname{A}\left(x\right) = \operatorname{A}\left(y\right)$ and $x \ne y$} \\
      \frac{P_{xy}}{\omega} & \mbox{if $\operatorname{A}\left(x\right) \ne \operatorname{A}\left(y\right)$,} \\
      -\sum\limits_{z \ne x} \frac{\partial P_{xz}}{\partial \omega} & \mbox{if $x = y$.}
      \end{cases}

Substitution model stationary state and derivatives
---------------------------------------------------------
The stationary state of the substitution model defined by :math:`P_{xy}` is

.. math::
  :label: px_M0

  p_{x} = \Phi_x

The derivatives of the stationary state with respect to :math:`\kappa` and :math:`\omega` are zero as these do not affect that state, so:

.. math::
   :label: dprx_dkappadomega_M0

   \frac{\partial p_{x}}{\partial \kappa} = \frac{\partial p_{x}}{\partial \omega} = 0

Calculation of induced equilibrium nucleotide frequencies
------------------------------------------------------

Following Pond *et al* 2010 (corrected), we calculate the codon frequencies :math:`\Phi_x` from the observed nucleotide frequencies, :math:`e^p_w`.

Let :math:`e^p_w` be the frequency of nucleotide :math:`w` at codon position :math:`p` in the alignment.

:math:`\Phi_x = \phi^1_{x_1} \times \phi^2_{x_2} \times \phi^3_{x_3}` where codon :math:`x` has nucleotide :math:`x_1` at codon position :math:`1`, nucleotide :math:`x_2` at codon position :math:`2`,
and nucleotide :math:`x_3` at codon position :math:`3`.
We can calculate :math:`\phi^p_w` by solving the following three equations:

.. math::
      :label: phi_pw

      e^1_w = \frac{\phi^1_w \times \left(1- \sum\limits_{wyz\epsilon X} \phi^2_y\times\phi^3_z\right)}{1-\pi_X}\\
      e^2_w = \frac{\phi^2_w \times \left(1- \sum\limits_{ywz\epsilon X} \phi^1_y\times\phi^3_z\right)}{1-\pi_X}\\
      e^3_w = \frac{\phi^3_w \times \left(1- \sum\limits_{yzw\epsilon X} \phi^1_y\times\phi^2_z\right)}{1-\pi_X}\\

where :math:`X = \{TAA, TAG, TGA\}` and :math:`\pi_X = \sum\limits_{xyz\epsilon X} \phi^1_x\times\phi^2_y\times\phi^3_z`.
