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
   \frac{P_{r,xy}}{\kappa} & \mbox{if $x$ is converted to $y$ by a transition of a nucleotide to $w$,} \\
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
   = \frac{P_{r,xy} \left(1 - \frac{F_{r,xy}}{\omega} \left(\pi_{r,\operatorname{A}\left(x\right)} / \pi_{r,\operatorname{A}\left(y\right)}\right)^{\beta}\right)}{\beta}
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


.. include:: weblinks.txt
