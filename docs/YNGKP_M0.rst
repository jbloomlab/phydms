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
