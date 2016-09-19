.. _implementation:

=======================================================
Numerical implementation of likelihood and gradients
=======================================================

.. contents::
   :depth: 2

*ExpCM* substitution model parameters and derivatives
---------------------------------------------------------
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

We define a variable transformation of the four nucleotide frequency parameters :math:`\phi_w` (three of which are unique). 
This transformation aids in numerical optimization. 
Specifically, we number the four nucleotides in alphabetical order so that :math:`w = 0` denotes *A*, :math:`w = 1` denotes *C*, :math:`w = 2` denotes *G*, and :math:`w = 3` denotes *T*. 
We then define the three free variables :math:`\eta_0`, :math:`\eta_1`, and :math:`\eta_2`, all of which are constrained to fall between zero and one.
For notational convenience in the formulas below, we also define :math:`\eta_3 = 0` -- note however that :math:`\eta_3` is not a free parameter, as it is always zero.
We define :math:`\phi_w` in terms of these :math:`\eta_i` variables by

.. math::
   :label: phi_from_eta

   \phi_w = \left(\prod_{i = 0}^{w - 1} \eta_i\right) \left(1 - \eta_w\right)

or conversely

.. math::
   :label: eta_from_phi

   \eta_w = 1 - \phi_w / \left(\prod_{i = 0}^{w - 1} \eta_i\right).

Note that setting :math:`\eta_w = \frac{3 - w}{4 - w}` makes all of the :math:`\phi_w` values equal to :math:`\frac{1}{4}`.

The derivatives are:

.. math::
   :label: dphi_deta

   \frac{\partial \phi_w}{\partial \eta_i} 
   & = & 
   \begin{cases}
   \left(\prod_{j = 0}^{i - 1} \eta_j\right)\left(\prod_{j = i + 1}^{w - 1}\eta_j\right) \left(1 - \eta_w\right) = \frac{\phi_w}{\eta_i} & \mbox{if $i < w$} \\
   -\prod_{j = 0}^{w - 1} \eta_j = \frac{\phi_w}{\eta_i - 1}& \mbox{if $i = w$} \\
   0 & \mbox{if $i > w$} \\
   \end{cases} \\
   & = & 
   \begin{cases}
   \frac{\phi_{w}}{\eta_i - \delta_{iw}} & \mbox{if $i \le w$,} \\
   0 & \mbox{otherwise,}
   \end{cases}

where :math:`\delta_{ij}` is the Kronecker-delta, equal to 1 if :math:`i = j` and 0 otherwise.


Given these definitions, the free parameters in and *ExpCM* model are :math:`\kappa`, :math:`\eta_0`, :math:`\eta_1`, :math:`\eta_2`, :math:`\beta`, and :math:`\omega`.

Here are the derivatives of :math:`P_{r,xy}` with respect to each of these parameters:

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
   :label: dPrxy_deta

   \frac{\partial P_{r,xy}}{\partial \eta_i}
   =
   \begin{cases}
   \frac{P_{r,xy}}{\phi_w} \frac{\partial \phi_w}{\partial \eta_i} = \frac{P_{r,xy}}{\eta_i - \delta_{iw}} & \mbox{if $x$ is converted to $y$ by a single-nucleotide mutation to $w \ge i$,} \\
   0 & \mbox{if $i > w$ or $x$ and $y$ differ by more than one nucleotide,} \\
   -\sum_z \frac{\partial P_{r,xz}}{\partial \eta_i} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_domega

   \frac{\partial P_{r,xy}}{\partial \omega} =
   \begin{cases}
   0 & \mbox{if $\operatorname{A}\left(x\right) = \operatorname{A}\left(y\right)$ and $x \ne y$} \\
   \frac{P_{r,xy}}{\omega} & \mbox{if $\operatorname{A}\left(x\right) \ne \operatorname{A}\left(y\right)$,} \\
   -\sum_z \frac{\partial P_{r,xy}}{\partial \omega} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_dbeta

   \frac{\partial P_{r,xy}}{\partial \beta} =
   \begin{cases}
   0 & \mbox{if $\operatorname{A}\left(x\right) = \operatorname{A}\left(y\right)$ and $x \ne y$}, \\
   0 & \mbox{if $\pi_{r,\operatorname{A}\left(x\right)} = \pi_{r,\operatorname{A}\left(y\right)}$ and $x \ne y$}, \\
   \frac{P_{r,xy}}{\beta} + P_{r,xy} 
   \frac{\left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta} \times \ln\left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)}{1 - \left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta}} 
   = \frac{P_{r,xy} \left(1 - F_{r,xy} \left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta}\right)}{\beta}
   & \mbox{if $\operatorname{A}\left(x\right) \ne \operatorname{A}\left(y\right)$,} \\
   -\sum_z \frac{\partial P_{r,xy}}{\partial \beta} & \mbox{if $x = y$.}
   \end{cases}

*ExpCM* stationary state and derivatives
------------------------------------------
The stationary state of the substitution model defined by :math:`P_{r,xy}` is

.. math::
   :label: prx

   p_{r,x} = \frac{q_x f_{r,x}}{\sum_z q_z f_{r,z}}

where

.. math::
   :label: frx

   f_{r,x} = \left(\pi_{r,\operatorname{A}\left(x\right)}\right)^{\beta}
   
and

.. math::
   :label: qx
   
   q_x = \phi_{x_0}\phi_{x_1}\phi_{x_2}

where :math:`x_0`, :math:`x_1`, and :math:`x_2` are the nucleotides at the first, second, and third positions of codon :math:`x`.

The derivatives of the stationary state with respect to :math:`\kappa` and :math:`\omega` are zero as these do not affect that state, so:

.. math::
   :label: dprx_dkappadomega

   \frac{\partial p_{r,x}}{\partial \kappa} = \frac{\partial p_{r,x}}{\partial \omega} = 0

The stationary state is sensitive to the value of :math:`\beta`, with derivative:

.. math::
   :label: dprx_dbeta

   \frac{\partial p_{r,x}}{\partial \beta} =
   \frac{p_{r,x}\left[\ln\left(\pi_{r,\operatorname{A}\left(x\right)}\right)\left(\sum_z f_{r,z} q_z\right) - \sum_z \ln\left(\pi_{r,\operatorname{A}\left(z\right)}\right) f_{r,z} q_z\right]}
   {\sum_z q_z f_{r,z}}

The stationary state is also sensitive to the values of :math:`\eta_0`, :math:`\eta_1`, and :math:`\eta_2`:

.. math::
   :label: dprx_detai

   \frac{\partial p_{r,x}}{\partial \eta_i} &=&
   \frac{f_{r,x} \frac{\partial q_x}{\partial \eta_i}\left(\sum_z q_z f_{r,z}\right) - f_{r,x} q_x \left(\sum_z f_{r,z} \frac{\partial q_z}{\partial \eta_i} \right)}
   {\left(\sum_z q_z f_{r,z}\right)^2} \\
   &=& \frac{\partial q_x}{\partial \eta_i}\frac{p_{r,x}}{q_x} - p_{r,x} \frac{\sum_z f_{r,z} \frac{\partial q_z}{\partial \eta_i}}{\sum_z q_z f_{r,z}} 

where the :math:`\frac{\partial q_x}{\partial \eta_i}` terms are:

.. math::
   :label: dqx_deta

   \frac{\partial q_x}{\partial \eta_i} & = & 
   \frac{\partial \phi_{x_0}}{\partial \eta_i} \phi_{x_1} \phi_{x_2} +
   \frac{\partial \phi_{x_1}}{\partial \eta_i} \phi_{x_0} \phi_{x_2} +
   \frac{\partial \phi_{x_2}}{\partial \eta_i} \phi_{x_0} \phi_{x_1} \\
   & = & 
   \sum_{j=0}^{2} \frac{\partial \phi_{x_j}}{\partial \eta_i} \prod_{k \ne j} \phi_{x_k} \\
   & = &
   q_x \sum_{j=0}^{2} \frac{1}{\phi_{x_j}} \frac{\partial \phi_{x_j}}{\partial \eta_i} \\
   & = &
   q_x \sum_{j=0}^{2} \frac{\operatorname{bool}\left(i \le x_j \right)}{\eta_i - \delta_{ix_j}}

where :math:`\operatorname{bool}\left(i \le j\right)` is 1 if :math:`i \le j` and 0 otherwise.

.. include:: weblinks.txt
