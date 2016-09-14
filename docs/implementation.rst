.. _implementation:

=======================================================
Numerical implementation of likelihood and gradients
=======================================================

.. contents::
   :depth: 2

Basic *ExpCM* substitution model
------------------------------------
We begin by considering the basic *ExpCM* substitution model defined in :ref:`ExpCM`. 
:math:`P_{r,xy}` gives the substitution rate from codon :math:`x` to codon :math:`y \ne x` at site :math:`r`, and is defined by

.. math::
   :label: Prxy

   P_{r,xy} = 
   \begin{cases}
   Q_{xy} \times F_{r,xy} & \mbox{if $x \ne y$,} \\
   -\sum_z F_{r,xz} Q_{xz} & \mbox{if $x = y$.}
   \end{cases}

where :math:`Q_{xy}` is the rate of mutation from codon :math:`x` to :math:`y` and is defined by

.. math::
   :label: Qxy

   Q_{xy} 
   &=& 
   \begin{cases} 
   \phi_w & \mbox{if $x$ is converted to $y$ by a single-nucleotide transversion to $w$,} \\
   \kappa \phi_w & \mbox{if $x$ is converted to $y$ by a single-nucleotide transition to $w$,} \\
   0 & \mbox{if $x$ and $y$ differ by more than one nucleotide,} \\
   \end{cases}

where :math:`\kappa` is the transition-transversion ratio and :math:`\phi_w` is the expected frequency of nucleotide :math:`w` in the absence of selection on amino-acid mutation (and so is subject to the constraint :math:`1 = \sum_w \phi_w`).

The "fixation probability" :math:`F_{r,xy}` of the mutation from :math:`x` to :math:`y` is

.. math::
   :label: Frxy

   F_{r,xy} 
   &=&
   \begin{cases}
   1 & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
   \omega & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\pi_{r,\mathcal{A}\left(x\right)} = \pi_{r,\mathcal{A}\left(y\right)}$} \\
   \omega \times \frac{-\beta \ln\left(\pi_{r,\mathcal{A}\left(x\right)} / \pi_{r,\mathcal{A}\left(y\right)}\right)}{1 - \left(\pi_{r,\mathcal{A}\left(x\right)} / \pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta}} & \mbox{otherwise.}
   \end{cases}

where :math:`\pi_{r,a}` is the preference of site :math:`r` for amino acid :math:`a`, :math:`\operatorname{A}\left(x\right)` is the amino acid encoded by codon :math:`x`, :math:`\beta` is the *stringency parameter*, and :math:`\omega` is a relative rate of nonsynonymous to synonymous mutations after accounting for the selection encapsulated by the preferences.

We define a variable transformation of the four the :math:`\phi_w` values (three of which are unique) to three transformed variables that are constrained to fall between zero and one. This transformation aids in numerical optimization.  Specifically, define:

.. math::
   :label: eta_from_phi
   
   \eta_0 &=& \phi_C + \phi_G \\
   \eta_1 &=& \frac{\phi_C}{\phi_C + \phi_G} \\
   \eta_2 &=&  \frac{\phi_A}{\phi_A + \phi_T}

or conversely:

.. math::
   :label: phi_from_eta

   \phi_A &=& \eta_2 \times \left(1 - \eta_0\right) \\
   \phi_C &=& \eta_1 \times \eta_0 \\
   \phi_G &=&  \left(1 - \eta_1\right) \times \eta_0 \\
   \phi_T &=&  \left(1 - \eta_2\right) \times \left(1 - \eta_0\right) \\

The free parameters in this model are :math:`\kappa`, :math:`\eta_0`, :math:`\eta_1`, :math:`\eta_2`, :math:`\beta`, and :math:`\omega`.

Here are the derivatives of :math:`P_{r,xy}` with respect to each of these parameters are:

.. math::
   :label: dPrxy_dkappa

   \frac{\partial P_{r,xy}}{\partial \kappa} 
   &=& 
   \begin{cases}
   \phi_w \times F_{r,xy} & \mbox{if $x$ is converted to $y$ by a transition of a nucleotide to $w$,} \\
   0 & \mbox{if $x$ and $y$ differ by something other than a single transition,} \\
   -\sum_z \frac{\partial P_{r,xz}}{\partial \kappa} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_deta0

   \frac{\partial P_{r,xy}}{\partial \eta_0} 
   &=&
   \begin{cases}
   - F_{r,xy} \eta_2 & \mbox{if $x$ is converted to $y$ by a transversion to A,} \\
   F_{r,xy} \eta_1 & \mbox{if $x$ is converted to $y$ by a transversion to C,} \\
   - F_{r,xy} \eta_1 & \mbox{if $x$ is converted to $y$ by a transversion to G,} \\
   F_{r,xy} \left(\eta_2 - 1\right) & \mbox{if $x$ is converted to $y$ by a transversion to T,} \\
   - \kappa F_{r,xy} \eta_2 & \mbox{if $x$ is converted to $y$ by a transition to A,} \\
   \kappa F_{r,xy} \eta_1 & \mbox{if $x$ is converted to $y$ by a transition to C,} \\
   - \kappa F_{r,xy} \eta_1 & \mbox{if $x$ is converted to $y$ by a transition to G,} \\
   \kappa F_{r,xy} \left(\eta_2 - 1\right) & \mbox{if $x$ is converted to $y$ by a transition to T,} \\
   0 & \mbox{otherwise if $x \ne y$},\\
   -\sum_z \frac{\partial P_{r,xz}}{\partial \eta_0} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_deta1

   \frac{\partial P_{r,xy}}{\partial \eta_1} 
   &=&
   \begin{cases}
   F_{r,xy} \eta_0 & \mbox{if $x$ is converted to $y$ by a transversion to C,} \\
   - F_{r,xy} \eta_0 & \mbox{if $x$ is converted to $y$ by a transversion to G,} \\
   \kappa F_{r,xy} \eta_0 & \mbox{if $x$ is converted to $y$ by a transition to C,} \\
   - \kappa F_{r,xy} \eta_0 & \mbox{if $x$ is converted to $y$ by a transition to G,} \\
   0 & \mbox{otherwise if $x \ne y$},\\
   -\sum_z \frac{\partial P_{r,xz}}{\partial \eta_1} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_deta2

   \frac{\partial P_{r,xy}}{\partial \eta_2} 
   &=&
   \begin{cases}
   F_{r,xy} \left(1 - \eta_0\right) & \mbox{if $x$ is converted to $y$ by a transversion to A,} \\
   F_{r,xy} \left(\eta_0 - 1\right) & \mbox{if $x$ is converted to $y$ by a transversion to T,} \\
   \kappa F_{r,xy} \left(1 - \eta_0\right) & \mbox{if $x$ is converted to $y$ by a transition to A,} \\
   \kappa F_{r,xy} \left(\eta_0 - 1\right) & \mbox{if $x$ is converted to $y$ by a transition to T,} \\
   0 & \mbox{otherwise if $x \ne y$},\\
   -\sum_z \frac{\partial P_{r,xz}}{\partial \eta_2} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_domega

   \frac{\partial P_{r,xy}}{\partial \omega} =
   \begin{cases}
   0 & \mbox{if $\operatorname{A}\left(x\right) = \operatorname{A}\left(y\right)$ and $x \ne y$} \\
   Q_{xy} & \mbox{if $\operatorname{A}\left(x\right) \ne \operatorname{A}\left(y\right)$ and $\pi_{r,\operatorname{A}\left(x\right)} = \pi_{r,\operatorname{A}\left(y\right)}$,} \\
   Q_{xy} \times \frac{-\beta \ln\left(\pi_{r,\mathcal{A}\left(x\right)} / \pi_{r,\mathcal{A}\left(y\right)}\right)}{1 - \left(\pi_{r,\mathcal{A}\left(x\right)} / \pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta}} & \mbox{if $\operatorname{A}\left(x\right) \ne \operatorname{A}\left(y\right)$ and $\pi_{r,\operatorname{A}\left(x\right)} \ne \pi_{r,\operatorname{A}\left(y\right)}$,} \\
   -\sum_z \frac{\partial P_{r,xy}}{\partial \omega} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_dbeta

   \frac{\partial P_{r,xy}}{\partial \beta} =
   \begin{cases}
   0 & \mbox{if $\operatorname{A}\left(x\right) = \operatorname{A}\left(y\right)$ and $x \ne y$}, \\
   0 & \mbox{if $\pi_{r,\operatorname{A}\left(x\right)} = \pi_{r,\operatorname{A}\left(y\right)}$ and $x \ne y$}, \\
   Q_{xy} \omega \left[\frac{-\ln\left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)}{1 - \left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta}} - \frac{\beta \times \left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta} \times \left[\ln\left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)\right]^2}{\left[1 - \left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta}\right]^2}\right] & \mbox{if $\operatorname{A}\left(x\right) \ne \operatorname{A}\left(y\right)$,} \\
   -\sum_z \frac{\partial P_{r,xy}}{\partial \beta} & \mbox{if $x = y$.}
   \end{cases}

.. include:: weblinks.txt
