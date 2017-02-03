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
   -\sum\limits_{z \ne x} F_{r,xz} Q_{xz} & \mbox{if $x = y$.}
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
   \frac{P_{r,xy}}{\kappa} & \mbox{if $x$ is converted to $y$ by a transition of a nucleotide to $w$,} \\
   0 & \mbox{if $x$ and $y$ differ by something other than a single transition,} \\
   -\sum\limits_{z \ne x} \frac{\partial P_{r,xz}}{\partial \kappa} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_deta

   \frac{\partial P_{r,xy}}{\partial \eta_i}
   =
   \begin{cases}
   \frac{P_{r,xy}}{\phi_w} \frac{\partial \phi_w}{\partial \eta_i} = \frac{P_{r,xy}}{\eta_i - \delta_{iw}} & \mbox{if $x$ is converted to $y$ by a single-nucleotide mutation to $w \ge i$,} \\
   0 & \mbox{if $i > w$ or $x$ and $y$ differ by more than one nucleotide,} \\
   -\sum\limits_{z \ne x} \frac{\partial P_{r,xz}}{\partial \eta_i} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_domega

   \frac{\partial P_{r,xy}}{\partial \omega} =
   \begin{cases}
   0 & \mbox{if $\operatorname{A}\left(x\right) = \operatorname{A}\left(y\right)$ and $x \ne y$} \\
   \frac{P_{r,xy}}{\omega} & \mbox{if $\operatorname{A}\left(x\right) \ne \operatorname{A}\left(y\right)$,} \\
   -\sum\limits_{z \ne x} \frac{\partial P_{r,xz}}{\partial \omega} & \mbox{if $x = y$.}
   \end{cases}

.. math::
   :label: dPrxy_dbeta

   \frac{\partial P_{r,xy}}{\partial \beta} =
   \begin{cases}
   0 & \mbox{if $\operatorname{A}\left(x\right) = \operatorname{A}\left(y\right)$ and $x \ne y$}, \\
   0 & \mbox{if $\pi_{r,\operatorname{A}\left(x\right)} = \pi_{r,\operatorname{A}\left(y\right)}$ and $x \ne y$}, \\
   \frac{P_{r,xy}}{\beta} + P_{r,xy} 
   \frac{\left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta} \times \ln\left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)}{1 - \left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta}} 
   & \mbox{if $\operatorname{A}\left(x\right) \ne \operatorname{A}\left(y\right)$,} \\
   -\sum\limits_{z \ne x} \frac{\partial P_{r,xz}}{\partial \beta} & \mbox{if $x = y$.}
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

   \frac{\partial p_{r,x}}{\partial \beta} &=&
   \frac{p_{r,x}\left[\ln\left(\pi_{r,\operatorname{A}\left(x\right)}\right)\left(\sum_z f_{r,z} q_z\right) - \sum_z \ln\left(\pi_{r,\operatorname{A}\left(z\right)}\right) f_{r,z} q_z\right]}
   {\sum_z q_z f_{r,z}} \\
   &=& p_{r,x} \left(\ln\left(\pi_{r,\operatorname{A}\left(x\right)}\right) - \frac{\sum_z \ln\left(\pi_{r,\operatorname{A}\left(z\right)}\right) f_{r,z} q_z}{\sum_z q_z f_{r,z}}\right) \\
   &=& p_{r,x} \left(\ln\left(\pi_{r,\operatorname{A}\left(x\right)}\right) - \sum_z \ln\left(\pi_{r,\operatorname{A}\left(z\right)}\right) p_{r,z}\right)

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

where :math:`\operatorname{bool}\left(i \le j\right)` is 1 if :math:`i \le j` and 0 otherwise, and so

.. math::
   :label: dprx_detai_2

   \frac{\partial p_{r,x}}{\partial \eta_i} & = &
   p_{r,x} 
   \left[\sum_{j=0}^{2} \frac{\operatorname{bool}\left(i \le x_j \right)}{\eta_i - \delta_{ix_j}} - \frac{\sum_z f_{r,z} q_z \sum_{j=0}^{2} \frac{\operatorname{bool}\left(i \le z_j \right)}{\eta_i - \delta_{iz_j}}}{\sum_z q_z f_{r,z}} \right] \\
   & = &
   p_{r,x} 
   \left[\sum_{j=0}^{2} \frac{\operatorname{bool}\left(i \le x_j \right)}{\eta_i - \delta_{ix_j}} - \frac{\sum_z p_{r,z} \sum_{j=0}^{2} \frac{\operatorname{bool}\left(i \le z_j \right)}{\eta_i - \delta_{iz_j}}}{\sum_z p_{r,z}} \right]
   \\ 


ExpCM with empirical nucleotide frequencies
----------------------------------------------
In the description above, the nucleotide frequencies :math:`\phi_w` are fit as three free parameters. 
Now let's consider the case where we instead calculate them empirically to give a stationary state that implies nucleotide frequencies that match those empirically observed in the alignment.
This should be beneficial in terms of optimization because it reduces the number of model parameters that need to be optimized.

Let :math:`g_w` be the empirical frequency of nucleotide :math:`w` taken over all sites and sequences in the alignment. 
Obviously, :math:`1 = \sum_w g_w`.
We want to empirically set :math:`\phi_w` to some value :math:`\hat{\phi}_w` such that when :math:`q_x = \hat{\phi}_{x_0} \hat{\phi}_{x_1} \hat{\phi}_{x_2}` then

.. math::
   :label: phihat

   g_w & = & 
   \frac{1}{L} \sum_r \sum_x \frac{1}{3} N_w\left(x\right) p_{r,x} \\ 
   & = & 
   \frac{1}{3L} \sum_r \frac{\sum_x N_w\left(x\right) f_{r,x} \prod_{k=0}^2 \hat{\phi}_{x_k}}{\sum_y f_{r,y} \prod_{k=0}^2 \hat{\phi}_{y_k}} \\

where :math:`N_w\left(x\right) = \sum_{k=0}^2 \delta_{x_k, w}` is the number of occurrence of nucleotide :math:`w` in codon :math:`x`, :math:`r` ranges over all codon sites in the gene, :math:`x` ranges over all codons, and :math:`k` ranges over the first three nucleotides.

There are three independent :math:`g_w` values and three independent :math:`\hat{\phi}_w` values (since :math:`1 = \sum_w g_w \sum_w = \hat{\phi}_w`), so we have three equations and three unknowns.
We could not solve the set of three equations analytically for the :math:`\hat{\phi}_w` values, so instead we use a non-linear equation solver to determine their values.

When using empirical nucleotide frequencies, we no longer need to calculate any derivatives with respect to :math:`\eta_i` as we no longer have the :math:`\eta_i` free parameters.

However, now the value of :math:`\phi_w = \hat{\phi}_w` depends on :math:`\beta` via the :math:`f_{r,x}` parameters in Equation :eq:`phihat`.
So we need new formulas for :math:`\frac{\partial p_{r,x}}{\partial \beta}` and :math:`\frac{\partial P_{r,xy}}{\partial \beta}` that accounts for this dependency.

Since we do not have an analytic expression for :math:`\hat{\phi}_w`, we cannot compute :math:`\frac{\partial \phi_w}{\partial \beta}` analytically.
But we can compute these derivatives numerically.
This is done using a finite-difference method.

We now update the formula in Equation :eq:`dPrxy_dbeta` for the case when `\phi_w` depends on `\beta`. In that case, we have:

.. math::
   :label: dQxy_dbeta_empirical_phi

   \frac{\partial Q_{xy}}{\partial \beta} = 
   \begin{cases}
   \frac{\partial \phi_w}{\partial \beta} & \mbox{if $x$ is converted to $y$ by a single-nucleotide transversion to $w$,} \\
   \kappa \frac{\partial \phi_w}{\partial \beta} & \mbox{if $x$ is converted to $y$ by a single-nucleotide transition to $w$,} \\
   0 & \mbox{if $x$ and $y$ differ by more than one nucleotide,} \\
   \end{cases}

and 

.. math::
   :label: dFrxy_dbeta_empirical_phi

   \frac{\partial F_{r,xy}}{\partial \beta}
   &=&
   \begin{cases}
   0 & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
   0 & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\pi_{r,\mathcal{A}\left(x\right)} = \pi_{r,\mathcal{A}\left(y\right)}$} \\
   \frac{F_{r,xy} \left(1 - \frac{F_{r,xy}}{\omega} \left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta}\right)}{\beta}
   & \mbox{otherwise,}
   \end{cases}

so for all :math:`x \ne y`, we have

.. math::
   :label: dPrxy_dbeta_empirical_phi

   \frac{\partial P_{r,xy}}{\partial \beta}
   & = &
   \frac{\partial \left(Q_{xy} \times F_{r,xy}\right)}{\partial \beta} \\
   & = &
   Q_{xy} \frac{\partial F_{r,xy}}{\partial \beta} + F_{r,xy} \frac{\partial Q_{xy}}{\partial \beta} \\
   & = &
   \left[\frac{\partial P_{r,xy}}{\partial \beta}\right]_{\mbox{free } \phi_w} + F_{r,xy} \frac{\partial Q_{xy}}{\partial \beta}.

where :math:`\left[\frac{\partial P_{r,xy}}{\partial \beta}\right]_{\mbox{free } \phi_w}` is the expression given by Equation :eq:`dPrxy_dbeta`.
When :math:`x = y`, we have :math:`\frac{\partial P_{r,xx}}{\partial \beta} = \sum\limits_{z \ne x} -\frac{\partial P_{r,xz}}{\partial \beta}`. 

We also must update the formula in Equation :eq:`dprx_dbeta` for the case where :math:`\phi_w` depends on :math:`\beta`. 
We have:

.. math::
   :label: dqx_dbeta_empirical_phi

   \frac{\partial q_x}{\partial \beta} 
   & = &
   \frac{\partial \left(\phi_{x_0} \phi_{x_1} \phi_{x_2}\right)}{\partial \beta} \\
   & = &
   \sum\limits_{j=0}^2 \frac{\partial \phi_{x_j}}{\partial \beta} \prod_{k \ne j} \phi_{x_k} \\
   & = &
   q_x \sum\limits_{j=0}^2 \frac{1}{\phi_{x_j}} \frac{\partial \phi_{x_j}}{\partial \beta}

and 

.. math::
   :label: dfrx_dbeta_empirical_phi

   \frac{\partial f_{r,x}}{\partial \beta}
   =
   f_{r,x} \left[\ln\left(\pi_{r,\operatorname{A}\left(x\right)}\right)\right].

So:

.. math::
   :label: dprx_dbeta_empirical_phi

   \frac{\partial p_{r,x}}{\partial \beta}
   & = &
   \frac{\partial}{\partial \beta} \left(\frac{q_x f_{r,x}}{\sum_z q_z f_{r,z}}\right) \\
   & = &
   \frac{\left(q_x \frac{\partial f_{r,x}}{\partial \beta} + f_{r,x} \frac{\partial q_x}{\partial \beta}\right) \sum_z q_z f_{r,z} - q_x f_{r,x} \sum_z \left(q_z \frac{\partial f_{r,z}}{\partial \beta} + f_{r,z} \frac{\partial q_z}{\partial \beta}\right)}{\left(\sum_z q_z f_{r,z}\right)^2} \\
   & = &
   \frac{\left(q_x f_{r,x} \left[\ln\left(\pi_{r,\operatorname{A}\left(x\right)}\right)\right] + f_{r,x} q_x \sum\limits_{j=0}^2 \frac{1}{\phi_{x_j}} \frac{\partial \phi_{x_j}}{\partial \beta}\right) \sum_z q_z f_{r,z} - q_x f_{r,x} \sum_z \left(q_z f_{r,z} \left[\ln\left(\pi_{r,\operatorname{A}\left(z\right)}\right)\right] + f_{r,z} q_z \sum\limits_{j=0}^2 \frac{1}{\phi_{z_j}} \frac{\partial \phi_{z_j}}{\partial \beta}\right)}{\left(\sum_z q_z f_{r,z}\right)^2} \\
   & = &
   p_{r,x} \frac{\left(\left[\ln\left(\pi_{r,\operatorname{A}\left(x\right)}\right)\right] + \sum\limits_{j=0}^2 \frac{1}{\phi_{x_j}} \frac{\partial \phi_{x_j}}{\partial \beta}\right) \sum_z q_z f_{r,z} - \sum_z \left(q_z f_{r,z} \left[\ln\left(\pi_{r,\operatorname{A}\left(z\right)}\right)\right] + f_{r,z} q_z \sum\limits_{j=0}^2 \frac{1}{\phi_{z_j}} \frac{\partial \phi_{z_j}}{\partial \beta}\right)}{\sum_z q_z f_{r,z}} \\
   & = &
   p_{r,x} \left[\left[\ln\left(\pi_{r,\operatorname{A}\left(x\right)}\right)\right] + \sum\limits_{j=0}^2 \frac{1}{\phi_{x_j}} \frac{\partial \phi_{x_j}}{\partial \beta} - \sum_z p_{r,z}\left(\left[\ln\left(\pi_{r,\operatorname{A}\left(z\right)}\right)\right] + \sum\limits_{j=0}^2 \frac{1}{\phi_{z_j}} \frac{\partial \phi_{z_j}}{\partial \beta}\right)\right] \\
   & = &
   \left[\frac{\partial p_{r,x}}{\partial \beta}\right]_{\mbox{free } \phi_w} +
   p_{r,x} \left[\sum\limits_{j=0}^2 \frac{1}{\phi_{x_j}} \frac{\partial \phi_{x_j}}{\partial \beta} - \sum_z p_{r,z}\sum\limits_{j=0}^2 \frac{1}{\phi_{z_j}} \frac{\partial \phi_{z_j}}{\partial \beta}\right]

where :math:`\left[\frac{\partial p_{r,x}}{\partial \beta}\right]_{\mbox{free } \phi_w}` is the expresssion given by Equation :eq:`dprx_dbeta`.

ExpCM with empirical nucleotide frequencies and diversifying pressure
---------------------------------------------------------------------
The :math:`\omega` value in the previous models is the gene-wide relative rate of nonsynonymous to synonymous mutations after accounting for the differing preferences among sites.  
In some cases, it might be possible to specify *a priori* expections for the diversifying pressure at each site. 
For instance, viruses benefit from amino-acid change in sites targeted by the immune system and, consequently, these sites have a higher rate of amino-acid substitution than expected given their level of inherent functional constraint. 
We can incorporate our expectations for diversifying pressure at specific sites into the selection terms :math:`F_{r,xy}`.  

Let :math:`\delta_{r}` be the pre-determined diversifying pressure for amino-acid change at site :math:`r` in the protein. A large positive value of :math:`\delta_r` corresponds to high pressure for amino-acid diversification, and negative value corresponds to expected pressure against amino-acid diversification beyond that captured in the amino-acid preferences.  We then replace :math:`\omega` in Equation :eq:`Frxy` with the expression :math:`\omega\times\left(1+\omega_{2}\times\delta_{r}\right)`, resulting in selection terms: 

.. math::
   :label: Frxy_divpressure

   F_{r,xy} = 
   \begin{cases}
   1 & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
   \omega\times\left(1+\omega_{2}\times\delta_{r}\right) & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\pi_{r,\mathcal{A}\left(x\right)} = \pi_{r,\mathcal{A}\left(y\right)}$} \\
   \omega\times\left(1+\omega_{2}\times\delta_{r}\right) \times \frac{\ln\left(\left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta}\right)}{1 - \left(\left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta}\right)} & \mbox{otherwise.}
   \end{cases}  
   
Whereas before :math:`\omega` reflected the elevation of non-synonymous substitution rate (averaged across the entire gene) beyond that expected given the amino-acid preferences, now :math:`\omega` reflects a gene-wide rate of elevated non-synonymous substitution after taking into account the expected sites of diversifying pressure (as represented by :math:`\delta_r`) weighted by :math:`\omega_{2}\times\delta_{r}`. These new selection terms in equation Equation :eq:`Frxy_divpressure` are identical the selection terms in Equation :eq:`Frxy` when :math:`\omega_{2} = 0`.

To ensure a positive value of :math:`\omega\times\left(1+\omega_{2}\times\delta_{r}\right)`, we constrain :math:`\omega>0`, :math:`-1<\omega_2<\infty`, and :math:`\lvert\max_r\delta_r\rvert\le1`.  

We have added one more parameter to Equation :eq:`Prxy`, :math:`\omega_2`, so we need to add a new derivative, :math:`\frac{\partial P_{r,xy}}{\partial \omega_2}` :eq:`dPrxy_domega2_divpressure`.  

.. math::
   :label: dPrxy_domega2_divpressure

   \frac{\partial P_{r,xy}}{\partial \omega_2} =
   \begin{cases}
   0 & \mbox{if $\operatorname{A}\left(x\right) = \operatorname{A}\left(y\right)$ and $x \ne y$} \\
   \omega \times \delta_r \times \frac{\ln\left(\left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta}\right)}{1 - \left(\left(\pi_{r,\mathcal{A}\left(x\right)}\right)^{\beta} / \left(\pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta}\right)} \times Q_{xy} & \mbox{if $\operatorname{A}\left(x\right) \ne \operatorname{A}\left(y\right)$,} \\
   -\sum\limits_{z \ne x} \frac{\partial P_{r,xy}}{\partial \omega_2} & \mbox{if $x = y$.}
   \end{cases}. 

Exponentials of the substitution matrix and derivatives
--------------------------------------------------------
The definitions above can be used to define a set of matrices :math:`\mathbf{P_r} = \left[P_{r,xy}\right]` that give the rate of transition from codon :math:`x` to :math:`y` at site :math:`r`. A key computation is to compute the probability of a transition in some amount of elapsed time :math:`\mu t`. These probabilities are given by 

.. math::
   :label: Mr
   
   \mathbf{M_r}\left(\mu t\right) = e^{\mu t\mathbf{P_r}}.

In this section, we deal with how to compute :math:`\mathbf{M_r}\left(\mu t\right)` and its derivatives.

Because :math:`\mathbf{P_r}` is reversible with stationary state given by the vector :math:`\mathbf{p_r} = \left[p_{r,x}\right]`, then as described by `Bryant, Galtier, and Poursat (2005) <http://www.maths.otago.ac.nz/~dbryant/Papers/04IHPLikelihood.pdf>`_, the matrix :math:`\left[\operatorname{diag}\left(\mathbf{p_r}\right)\right]^{\frac{1}{2}} \mathbf{P_r} \left[\operatorname{diag}\left(\mathbf{p_r}\right)\right]^{\frac{-1}{2}}` is symmetric.

We can use a numerical routine to compute the eigenvalues and orthonormal eigenvectors. 
Let :math:`\mathbf{D_r}` be a diagonal matrix with elements equal to the eigenvalues, let :math:`\mathbf{B_r}` be the matrix whose columns are the right orthonormal eigenvectors (in the same order as the eigenvalues), and note that :math:`\mathbf{B_r}^{-1} = \mathbf{B_r}^T`.
Then we have :math:`\left[\operatorname{diag}\left(\mathbf{p_r}\right)\right]^{\frac{1}{2}} \mathbf{P_r} \left[\operatorname{diag}\left(\mathbf{p_r}\right)\right]^{\frac{-1}{2}} = \mathbf{B_r} \mathbf{D_r} \mathbf{B_r}^T` or equivalently

.. math::

   \mathbf{P_r} = \mathbf{A_r} \mathbf{D_r} \mathbf{A_r}^{-1}

where

.. math::

   \mathbf{A_r} = \left[\operatorname{diag}\left(\mathbf{p_r}\right)\right]^{\frac{-1}{2}} \mathbf{B_r}

and 

.. math::

   \mathbf{A_r}^{-1} = \mathbf{B_r}^T \left[\operatorname{diag}\left(\mathbf{p_r}\right)\right]^{\frac{1}{2}}.

The matrix exponentials are then easily calculated as

.. math::

   \mathbf{M_r}\left(\mu t\right) = e^{\mu t\mathbf{P_r}} = \mathbf{A_r} e^{\mu t \mathbf{D_r}} \mathbf{A_r}^{-1}

We also want to calculate the derivatives of :math:`\mathbf{M_r}\left(\mu t\right)` with respect to the other parameters on which :math:`P_{r,xy}` depends (e.g., :math:`\beta`, :math:`\eta_i`, :math:`\kappa`, and :math:`\omega`). 

According to `Kalbeisch and Lawless (1985) <https://www.jstor.org/stable/2288545?seq=1#page_scan_tab_contents>`_ (see also `Kenney and Gu (2012) <http://www.mathstat.dal.ca/~hgu/hession.pdf>`_ and `Bloom et al (2011) <http://dx.plos.org/10.1371/journal.pone.0022201>`_), the derivative with respect to some parameter :math:`z` is given by

.. math::

   \frac{\partial \mathbf{M_r}\left(\mu t\right)}{\partial z} = \mathbf{A_r} \mathbf{V_{r,z}} \mathbf{A_r}^{-1}

where the elements of :math:`\mathbf{V_{r,z}}` are

.. math::

   V_{xy}^{r,z} = 
   \begin{cases}
   B_{xy}^{r,z} \frac{\exp\left(\mu t D^r_{xx}\right) - \exp\left(\mu t D^r_{yy}\right)}{D^r_{xx} - D^r_{yy}} & \mbox{if $x \ne y$ and $D_{xx}^r \ne D_{yy}^r$}, \\
   B_{xy}^{r,z} \mu t \exp\left(\mu t D^r_{xx}\right) & \mbox{if $x \ne y$ and $D_{xx}^r = D_{yy}^r$}, \\
   B_{xx}^{r,z} \mu t \exp\left(\mu t D_{xx}^r\right)& \mbox{if $x = y$}, 
   \end{cases}

where :math:`D^r_{xx}` and :math:`D^r_{yy}` are the diagonal elements of :math:`\mathbf{D_r}`, and :math:`B_{xy}^{r,z}` are the elements of the matrix :math:`\mathbf{B_{r,z}}` defined by

.. math::
 
   \mathbf{B_{r,z}} = \mathbf{A_r}^{-1} \frac{\partial \mathbf{P_r}}{\partial z} \mathbf{A_r}.

Scaling the branch lengths with a mutation rate
-------------------------------------------------
The aforementioned section defines the substitution probabilities in terms of :math:`\mu t` (e.g., Eq. :eq:`Mr`). 
Here :math:`\mu` is a substitution rate, and :math:`t` is the branch length. 
If we are freely optimizing all branch lengths, then we just set :math:`\mu = 1` so that :math:`\mu t = t`, and then :math:`\mu` effectively drops out. 
However, if we have fixed the branch lengths are **not** optimizing them, then we might want to include a parameter :math:`\mu` that effectively re-scales all the fixed branch lengths by a constant. 
In this case, :math:`\mu` also becomes a free parameter of the model, and we want to compute the derivative of :math:`\mathbf{M_r}\left(\mu t\right)` with respect to :math:`\mu`. This is straightforward:

.. math::

   \frac{\partial \mathbf{M_r}\left(\mu t\right)}{\partial \mu} = t \mathbf{P_r} e^{\mu t \mathbf{P_r}}.


Calculating the likelihood and derivatives on a tree
------------------------------------------------------
Above we describe computing the transition probabilities as a function of branch length.
Here we consider how to use those computations to compute the actual likelihoods on a tree.

.. figure:: implementation_tree_example.jpg
   :align: center
   :alt: implementation_tree_example.jpg
   :width: 40%

   The tree used in the example calculation below.

We begin by computing the likelihood of the alignment at a specific site.
Let :math:`\mathcal{S}_r` denote the set of aligned codons at site :math:`r`, let :math:`\mathcal{T}` by the phylogenetic tree with branch lengths specified, and let :math:`\mathbf{P_r}` be the transition matrix at site :math:`r` defined above.
Then the likelihood at site :math:`r` is :math:`\Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right)`.
For the example tree above, we can use the pruning algorithm (`Felsenstein, J Mol Evol, 1981`_) to write

.. math::

   \Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right) 
   =
   \sum_y p_{r,y} M_{r,yGAA}\left(t_3\right) \left[\sum_x M_{r,yx}\left(t_4\right) M_{r,xCAA}\left(t_1\right) M_{r,xCAG}\left(t_2\right)\right].

Let :math:`n` denote a node on a tree, let :math:`t_n` denote the length of the branch leading to node :math:`n`, and let :math:`\mathcal{d}_1\left(n\right)` and :math:`\mathcal{d}_1\left(n\right)` denote the right and left descendents of node :math:`n` for all non-terminal nodes. 
Then define the *partial conditional likelihood* of the subtree rooted at :math:`n` as:

.. math::

   L_{r,n}\left(x\right) =
   \begin{cases}
   \delta_{x\mathcal{S}_{r,n}} & \mbox{if $n$ is a tip node with codon $\mathcal{S}_{r,n}$ at site $r$,} \\
   1 & \mbox{if $n$ is a tip node with a gap at site $r$,} \\
   \left[\sum_y M_{r,xy}\left(t_{\mathcal{d}_1\left(n\right)}\right) L_{r, \mathcal{d}_1\left(n\right)}\left(y\right)\right] \left[\sum_y M_{r,xy}\left(t_{\mathcal{d}_2\left(n\right)}\right) L_{r, \mathcal{d}_2\left(n\right)}\left(y\right)\right] & \mbox{otherwise.}
   \end{cases}

where :math:`\delta_{xy}` is the `Kronecker delta`_.
So for instance in the example tree above, :math:`L_{r,n_4}\left(x\right) = M_{r,xCAA}\left(t_1\right) M_{r,xCAG}\left(t_2\right)`, and :math:`L_{r,n_5}\left(y\right) = M_{r,yGAA} \sum_x M_{r,yx}\left(t_4\right) L_{r,n_4}\left(x\right)`.

Using this definition, we have

.. math::

   \Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right) 
   = \sum_x p_{r,x} L_{r,n_{\rm{root}}}\left(x\right)

where :math:`n_{\rm{root}}` is the root node of tree :math:`\mathcal{T}`; :math:`n_{\rm{root}} = n_5` in the example tree above.

In practice, we usually work with the log likelihoods (always using natural logarithms).
The total likelihood is the sum of the log likelihoods for each site:

.. math::

   \ln \left[ \Pr\left(\mathcal{S} \mid \mathcal{T}, \left\{\mathbf{P_r}\right\}\right) \right] = \sum_r \ln \left[\Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right) \right].

We next consider how to compute the derivatives with respect to some model parameter.
Let :math:`\alpha` denote the model parameter in question, and assume that we have already determined :math:`\frac{M_{r,xy}\left(t\right)}{\partial \alpha}`.
By the chain rule, we have

.. math::

   \frac{\partial L_{r,n}\left(x\right)}{\partial \alpha} =
   \begin{cases}
   0 & \mbox{if $n$ is a tip node,}, \\ 
   \\
   \left[\sum_y \frac{\partial M_{r,xy}\left(t_{\mathcal{d}_1\left(n\right)}\right)}{\partial \alpha} L_{r, \mathcal{d}_1\left(n\right)}\left(y\right) + M_{r,xy}\left(t_{\mathcal{d}_1\left(n\right)}\right) \frac{\partial L_{r, \mathcal{d}_1\left(n\right)}\left(y\right)}{\partial \alpha}\right] 
   \left[\sum_y M_{r,xy}\left(t_{\mathcal{d}_2\left(n\right)}\right) L_{r, \mathcal{d}_2\left(n\right)}\left(y\right)\right] & \mbox{otherwise.} \\
   + \left[\sum_y M_{r,xy}\left(t_{\mathcal{d}_1\left(n\right)}\right) L_{r, \mathcal{d}_1\left(n\right)}\left(y\right)\right] 
   \left[\sum_y \frac{\partial M_{r,xy}\left(t_{\mathcal{d}_2\left(n\right)}\right)}{\partial \alpha} L_{r, \mathcal{d}_2\left(n\right)}\left(y\right) + M_{r,xy}\left(t_{\mathcal{d}_2\left(n\right)}\right) \frac{\partial L_{r, \mathcal{d}_2\left(n\right)}\left(y\right)}{\partial \alpha}\right] & \\
   \end{cases}

The derivative of the likelihood at the site is then

.. math::

   \frac{\partial \Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right)}{\partial \alpha} 
   = \sum_x \left(\frac{\partial p_{r,x}}{\partial \alpha} L_{r,n_{\rm{root}}}\left(x\right) + p_{r,x} \frac{\partial L_{r,n_{\rm{root}}}\left(x\right)}{\partial \alpha}\right)

and the derivative of the log likelihood at the site is

.. math::

   \frac{\partial \ln\left[\Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right)\right]}{\partial \alpha} 
   = \frac{\sum_x \left(\frac{\partial p_{r,x}}{\partial \alpha} L_{r,n_{\rm{root}}}\left(x\right) + p_{r,x} \frac{\partial L_{r,n_{\rm{root}}}\left(x\right)}{\partial \alpha}\right)}
   {\Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right)}.

The derivative of the overall log likelihood is

.. math::

   \frac{\partial \ln \left[ \Pr\left(\mathcal{S} \mid \mathcal{T}, \left\{\mathbf{P_r}\right\}\right) \right]}{\partial \alpha} 
   = \sum_r \frac{\partial \ln \left[\Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right) \right]}{\partial \alpha}.

Scaling to avoid numerical underflow
--------------------------------------
For larger trees, there can be numerical underflow due to multiplication of lots of small numbers when computing the likelihoods. This issue, and how it can be solved by re-scaling the likelihoods during the calculation, is discussed on page 426 of `Yang, J Mol Evol, 51:423-432`_.

Let :math:`L_{r,n}\left(x\right)` be the partial conditional likelihood at node :math:`n` of codon :math:`x` at site :math:`r` as defined above. These partial conditional likelihoods can get very small as we move up the tree towards the root, as they are recursively defined as the products of very small numbers. For the scaling to avoid underflow, we define the scaled partial condition likelilhood as

.. math::

   \tilde{L}_{r,n}\left(x\right) = \frac{L_{r,n}\left(x\right)}{U_{r,n} \times \prod\limits_{k < n} U_{r,k}}

where we use :math:`k < n` to indicate all nodes :math:`k` that are descendants of :math:`n`, and where

.. math::

   U_{r,n} = 
   \begin{cases}
   1 & \mbox{if $n$ is divisible by $K$,} \\
   \max_x\left[L_{r,n}\left(x\right) \times \prod\limits_{k < n} U_{r,k}\right] & \mbox{otherwise}
   \end{cases}

where :math:`K` is the frequency with which we re-scale the likelihoods. A reasonable value of :math:`K` might be 5 or 10. Effectively, this means that every :math:`K` nodes we are re-scaling so that the largest partial conditional likelihood is one.

With this re-scaling, the total likelihood at site :math:`r` is then

.. math::

   \Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right)
   = \left(\sum\limits_x p_{r,x} \tilde{L}_{r,n_{\rm{root}}}\left(x\right)\right) \times
   \left(\prod\limits_n U_{r,n}\right)

and the total log likelihood at site :math:`r` is

.. math::

   \ln\left[\Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right)\right]
   = \ln\left(\sum\limits_x p_{r,x} \tilde{L}_{r,n_{\rm{root}}}\left(x\right)\right) +
   \sum\limits_n \ln\left(U_{r,n}\right).

The derivative is then

.. math::

   \frac{\partial \ln\left[\Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right)\right]}{\partial \alpha} &=&
   \frac{\frac{\partial}{\partial \alpha}\left[ \left(\sum\limits_x p_{r,x} \tilde{L}_{r,n_{\rm{root}}}\left(x\right)\right) \times \left(\prod\limits_n U_{r,n}\right) \right]}{\left(\sum\limits_x p_{r,x} \tilde{L}_{r,n_{\rm{root}}}\left(x\right)\right) \times \left(\prod\limits_n U_{r,n}\right)} \\
   &=&
   \frac{\left(\sum\limits_x \left[\frac{\partial p_{r,x}}{\partial \alpha} 
   \tilde{L}_{r,n_{\rm{root}}}\left(x\right) + p_{r,x}\frac{\partial \tilde{L}_{r,n_{\rm{root}}}\left(x\right)}{\partial \alpha}\right]\right) \times \left(\prod\limits_n U_{r,n}\right) + \left(\sum\limits_x p_{r,x} \tilde{L}_{r,n_{\rm{root}}}\left(x\right)\right) \times \frac{\partial \left(\prod\limits_n U_{r,n}\right)}{\partial \alpha}}
   {\left(\sum\limits_x p_{r,x} \tilde{L}_{r,n_{\rm{root}}}\left(x\right)\right) \times \left(\prod\limits_n U_{r,n}\right)} \\
   &=& 
   \frac{\sum_x \left(\frac{\partial p_{r,x}}{\partial \alpha} \tilde{L}_{r,n_{\rm{root}}}\left(x\right) + p_{r,x} \frac{\partial \tilde{L}_{r,n_{\rm{root}}}\left(x\right)}{\partial \alpha}\right)}{\sum\limits_x p_{r,x} \tilde{L}_{r,n_{\rm{root}}}\left(x\right)} + \frac{\frac{\partial \left(\prod\limits_n U_{r,n}\right)}{\partial \alpha}}{\prod\limits_n U_{r,n}}.

For reasons that are not immediately obvious to me but are clearly verified by numerical testing, this last term of :math:`\frac{\frac{\partial \left(\prod\limits_n U_{r,n}\right)}{\partial \alpha}}{\prod\limits_n U_{r,n}}` is zero, and so

.. math::

   \frac{\partial \ln\left[\Pr\left(\mathcal{S}_r \mid \mathcal{T}, \mathbf{P_r}\right)\right]}{\partial \alpha} =
   \frac{\sum_x \left(\frac{\partial p_{r,x}}{\partial \alpha} \tilde{L}_{r,n_{\rm{root}}}\left(x\right) + p_{r,x} \frac{\partial \tilde{L}_{r,n_{\rm{root}}}\left(x\right)}{\partial \alpha}\right)}{\sum\limits_x p_{r,x} \tilde{L}_{r,n_{\rm{root}}}\left(x\right)}.

In practice, we work with the :math:`\tilde{L}_{r,n}\left(x\right)` values, and then apply the correction of adding :math:`\sum_n \ln\left(U_r,n\right)` at the end.


Derivatives with respect to branch lengths
--------------------------------------------
The section above describes how to compute the derivatives with respect to paramters (e.g., model parameters) that affect all parts of the tree. 
In many cases, we may want to optimize individual branch lengths rather than the overall mutation rate :math:`\mu`.
In this case, we need to compute the derivatives with respect to the branch lengths.
This is somewhat simpler for each individual branch length, since each individual branch length only affects part of the tree.

Specifically, for each internal node :math:`n`,

.. math::

   \frac{\partial L_{r,n}\left(x\right)}{\partial t_{\mathcal{d}_1\left(n\right)}} =

Units of tree branch lengths
------------------------------
When we optimize with the :math:`P_{r,xy}` substitution matrices described above, the resulting branch lengths are **not** in units of substitutions per site. Therefore, for tree input / output, we re-scale the branch lengths so that they are in units of substitution per site.

In a single unit of time, the probability that if site :math:`r` is initially :math:`x`, then it will undergo a substitution to some other codon :math:`y` is :math:`\sum_{y \ne x} P_{r,xy} = -P_{r,xx}`. Since the equilibrium probability that site :math:`r` is :math:`x` is :math:`p_{r,x}`, then the probability that site :math:`r` undergoes a substitution in a unit of time is :math:`-\mu \sum_x p_{r,x} P_{r,xx}`. So averaging over all :math:`L` sites, the probability that the average site will undergo a substitution in a unit of time is :math:`-\frac{\mu}{L} \sum_{r=1}^L \sum_x p_{r,x} P_{r,xx}`.

Therefore, if we optimize the branch lengths :math:`t_b` and the model parameters in :math:`P_{r,xy}`, and then at the end re-scale the branch lengths to :math:`t_b' = t_b \times \frac{-\mu}{L} \sum_{r=1}^L \sum_x p_{r,x} P_{r,xx}` then the re-scaled branch lengths :math:`t_b` are in units of substitutions per sites. Therefore, for input and output to ``phydms``, we assume that input branch lengths are already in units of substitutions per site, and scale them from :math:`t_b'` to :math:`t_b`. Optimization is performed on :math:`t_b`, and then for output we re-scale the optimized branch lengths from :math:`t_b` to :math:`t_b'`.

.. include:: weblinks.txt
 


   
