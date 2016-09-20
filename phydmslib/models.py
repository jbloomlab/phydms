"""Substitution models.

For all models, nucleotides, amino acids, and codons are indexed
by integers 0, 1, ... using the indexing schemes defined
in `phydmslib.constants`.
"""


import scipy
from phydmslib.constants import *


class ExpCM:
    """A set of ExpCM for a gene.

    See the `__init__` method for how to initialize an `ExpCM`.

    Attributes should **only** be updated via the `updateParams`
    method, do **not** set attributes directly.

    The attributes are listed below. Note that only the first
    few represent independent parameters of the model. The rest
    are dependent attributes that are calculated from these
    independent parameters.

    Attributes: 
        `nsites` (int)
            Number of codon sites in gene.
        `pi` (`numpy.ndarray` of floats, shape `(nsites, N_AA)`)
            `pi[r][a]` is preference of site `r` for amino-acid `a`.
        `kappa` (float > 0)
            Transition-transversion ratio.
        `omega` (float > 0)
            Nonsynonymous to synonymous substitution ratio.
        `beta` (float > 0)
            Stringency parameter re-scaling amino-acid preferences.
        `phi` (`numpy.ndarray` of floats, length `N_NT`)
            Nucleotide frequency params summing to one, can be
            computed from `eta`.
        `eta` (`numpy.ndarray` of floats, length `N_NT - 1`)
            Transformation of the nucleotide frequency params in `phi`;
            all entries are > 0 and < 1.
        `Prxy` (`numpy.ndarray` of floats, shape `(nsites, N_CODON, N_CODON)`
            `Prxy[r][x][y]` is substitution rate from codon `x` to `y` at `r`.
            Diagonal elements make rows sum to zero for each `Prxy[r]`.
        `prx` (`numpy.ndarray` of floats, shape `(nsites, N_CODON)`
            `prx[r][x]` is stationary state of `Prxy` for codon `x` at `r`.
        `Qxy` (`numpy.ndarray` of floats, shape `(N_CODON, N_CODON)`
            `Qxy[x][y]` is mutation rate from `x` to `y`, diagonal undefined.
        `qx` (`numpy.ndarray` of floats, length `N_CODON`
            `qx[x]` is stationary state of `Qxy` for codon `x`.
        `Frxy` (`numpy.ndarray` of floats, shape `(nsites, N_CODON, N_CODON)`
            `Frxy[r][x][y]` fixation prob from `x` to `y`, diagonal 
            undefined for each `Frxy[r]`.
        `frx` (`numpy.ndarray` of floats, shape `(nsites, N_CODON)`
            `frx[r][x]` is stationary state of `Frxy` for codon `x` at `r`.
        `pi_codon` (`numpy.ndarray` of floats, shape `(nsites, N_CODON)`)
            `pi_codon[r][x]` is preference of site `r` for amino acid
            encoded by codon `x`.
        `piAx_piAy` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
            `piAx_piAy[r][x][y]` is ratio of preference for amino acid
            encoded by codon `x` divided by pref for that encoded by `y`.
        `piAx_piAy_beta` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
            Equal to `piAx_piAy` raised to the power of `beta`.
    """

    def __init__(self, prefs, kappa=1.0, omega=1.0, beta=1.0,
            phi=scipy.ones(N_NT) / N_NT):
        """Initialize an `ExpCM` object.

        Args: 
            `prefs` (list)
                List of dicts giving amino-acid preferences for
                each site. Each dict keyed by amino acid letter
                codes, value is pref > 0 and < 1. Must sum to 1 
                at each site.
            `kappa`, `omega`, `beta`, and `phi`
                Meanings described in main class doc string.
        """
        self.nsites = len(prefs)
        assert self.nsites > 0, "No preferences specified"

        # put prefs in pi, adjust to > ALMOST_ZERO
        self.pi = scipy.ndarray((self.nsites, N_AA), dtype='float')
        assert (isinstance(prefs, list) and 
                all([isinstance(x, dict) for x in prefs])),\
                "prefs is not a list of dicts"
        for r in range(self.nsites):
            assert set(prefs[r].keys()) == set(AA_TO_INDEX.keys()),\
                    "prefs not keyed by amino acids for site {0}".format(r)
            assert abs(1 - sum(prefs[r].values())) <= ALMOST_ZERO,\
                    "prefs don't sum to one for site {0}".format(r)
            assert all([0 < pix < 1 for pix in prefs[r].values()]),\
                    "prefs aren't all > 0 and < 1 for site {0}".format(r)
            for (a, aa) in INDEX_TO_AA.items():
                self.pi[r][a] = prefs[r][aa]
            self.pi[self.pi < ALMOST_ZERO] = ALMOST_ZERO
            self.pi[r] /= self.pi[r].sum() # renormalize to sum to one 

        # set up attributes defined solely in terms of preferences
        self.pi_codon = scipy.full((self.nsites, N_CODON), -1, dtype='float')
        self._update_pi_codon()
        self.piAx_piAy = scipy.full((self.nsites, N_CODON, N_CODON), -1, 
                dtype='float')
        self._update_piAx_piAy()

        # construct eta from phi after ensuring > ALMOST_ZERO
        assert phi.dtype == 'float', 'phi not array of floats'
        assert len(phi) == N_NT, "not the right number of entries in phi"
        assert abs(1 - phi.sum()) <= ALMOST_ZERO, "phi doesn't sum to 1"
        assert (phi > 0).all() and (phi < 1).all(),\
                "phi entries not > 0 and < 1"
        self.phi = phi.copy()
        self.phi[self.phi < ALMOST_ZERO] = ALMOST_ZERO
        self.phi /= self.phi.sum()
        eta = scipy.ndarray(N_NT - 1, dtype='float')
        etaprod = 1.0
        for w in range(N_NT - 1):
            eta[w] = 1.0 - self.phi[w] / etaprod
            etaprod *= eta[w]

        # set attributes to calling params
        self.kappa = kappa
        self.omega = omega
        self.beta = beta
        self.eta = eta

        # define other params, initialized appropriately
        self.piAx_piAy_beta = scipy.zeros((self.nsites, N_CODON, N_CODON), 
                dtype='float')
        self.Prxy = scipy.zeros((self.nsites, N_CODON, N_CODON), dtype='float')
        self.prx = scipy.zeros((self.nsites, N_CODON), dtype='float')
        self.Qxy = scipy.zeros((N_CODON, N_CODON), dtype='float')
        self.qx = scipy.zeros(N_CODON, dtype='float')
        self.Frxy = scipy.ones((self.nsites, N_CODON, N_CODON), dtype='float')
        self.frx = scipy.zeros((self.nsites, N_CODON), dtype='float')
        
        self.updateParams(update_all=True)

    def updateParams(self, kappa=None, omega=None, beta=None, eta=None,
            update_all=False):
        """Update model params.

        This method is the **only** acceptable way to update `ExpCM`
        model parameters. Any arguments that are not set to the
        default value of `None` are updated to the new value. This
        method then updates any other dependent `ExpCM` attributes.

        Args:
            `kappa`, `omega`, `beta`, and `eta`
                Attribute meanings described in main class doc string.
            `update_all` (bool)
                If `True`, update all dependent attributes using
                current values of model parameters.

        """
        if update_all:
            changed = set(['kappa', 'eta', 'beta', 'omega'])
        else:
            changed = set([]) # contains string names of changed params

        if (eta is not None) and (eta != self.eta).any():
            assert (isinstance(eta, scipy.ndarray) and eta.shape == (3,) 
                    and eta.dtype == float), "eta not array of 3 floats"
            assert (eta > 0).all() and (eta < 1).all(),\
                    "eta must be > 0 and < 1"
            self.eta = eta
            changed.add('eta')

        if kappa != None and kappa != self.kappa:
            assert isinstance(kappa, float) and kappa > 0,\
                    "kappa must be float > 0"
            self.kappa = kappa
            changed.add('kappa')

        if omega != None and omega != self.omega:
            assert isinstance(omega, float) and omega > 0,\
                    "omega must be float > 0"
            self.omega = omega
            changed.add('omega')

        if beta != None and beta != self.beta:
            assert isinstance(beta, float) and beta > 0,\
                    "beta must be float > 0"
            self.beta = beta
            changed.add('beta')

        if 'eta' in changed:
            self._update_phi()
            self._update_Qxy()
            self._update_qx()
        elif 'kappa' in changed:
            self._update_Qxy()

        if 'beta' in changed:
            self._update_piAx_piAy_beta()
            self._update_Frxy()
            self._update_frx()
        elif 'omega' in changed:
            self._update_Frxy()

        if ('beta' in changed) or ('eta' in changed):
            self._update_Prxy()
            self._update_prx()
        elif changed:
            self._update_Prxy()

    def _update_phi(self):
        """Update `phi` using current `eta`."""
        etaprod = 1.0
        for w in range(N_NT - 1):
            self.phi[w] = etaprod * (1 - self.eta[w])
            etaprod *= self.eta[w]
        self.phi[N_NT - 1] = etaprod

    def _update_Qxy(self):
        """Update `Qxy` using current `kappa` and `phi`."""
        for w in range(N_NT):
            scipy.copyto(self.Qxy, self.phi[w], where=CODON_NT_MUT[w])
        self.Qxy[CODON_TRANSITION] *= self.kappa

    def _update_qx(self):
        """Update `qx` using current `phi`."""
        self.qx.fill(1.0)
        for j in range(3):
            for w in range(N_NT):
                self.qx[CODON_NT[j][w]] *= self.phi[w]

    def _update_pi_codon(self):
        """Update `pi_codon` using current `pi`."""
        for r in range(self.nsites):
            for a in range(N_AA):
                scipy.copyto(self.pi_codon[r], self.pi[r][a], 
                where=(CODON_TO_AA == a))

    def _update_piAx_piAy(self):
        """Update `piAx_piAy` from `pi_codon`."""
        with scipy.errstate(divide='raise', under='raise', over='raise', 
                invalid='raise'):
            for r in range(self.nsites):
                pim = scipy.tile(self.pi_codon[r], (N_CODON, 1)) # [x][y] is piAy
                scipy.copyto(self.piAx_piAy[r], pim.transpose() / pim)

    def _update_piAx_piAy_beta(self):
        """Update `piAx_piAy_beta` from `piAx_piAy` and `beta`."""
        with scipy.errstate(divide='raise', under='raise', over='raise', 
                invalid='raise'):
            scipy.copyto(self.piAx_piAy_beta, self.piAx_piAy**self.beta)

    def _update_Frxy(self):
        """Update `Frxy` from `piAx_piAy_beta`, `omega`, and `beta`."""
        with scipy.errstate(divide='raise', under='raise', over='raise', 
                invalid='ignore'):
            scipy.copyto(self.Frxy, -self.omega * scipy.log(
                    self.piAx_piAy_beta) / (1 - self.piAx_piAy_beta), 
                    where=CODON_NONSYN)
        scipy.copyto(self.Frxy, self.omega, where=(1 == self.piAx_piAy_beta))
        assert not scipy.isnan(self.Frxy).any(), "NaN values in Frxy"
        assert not scipy.isinf(self.Frxy).any(), "Infinite values in Frxy"

    def _update_frx(self):
        """Update `frx` using current `pi_codon` and `beta`."""
        scipy.copyto(self.frx, self.pi_codon**self.beta)

    def _update_Prxy(self):
        """Update `Prxy` using current `Frxy` and `Qxy`."""
        scipy.copyto(self.Prxy, self.Frxy * self.Qxy)
        diag_indices = scipy.diag_indices(N_CODON)
        for r in range(self.nsites):
            scipy.fill_diagonal(self.Prxy[r], 0)
            self.Prxy[r][diag_indices] = -scipy.sum(self.Prxy[r], axis=1)

    def _update_prx(self):
        """Update `prx` using current `frx` and `qx`."""
        self.prx = self.frx * self.qx
        with scipy.errstate(divide='raise', under='raise', over='raise',
                invalid='raise'):
            for r in range(self.nsites):
                self.prx[r] /= self.prx[r].sum()



if __name__ == '__main__':
    import doctest
    doctest.testmod()
