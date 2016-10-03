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
            Nucleotide frequency params summing to one, is
            computed from `eta`.
        `eta` (`numpy.ndarray` of floats, length `N_NT - 1`)
            Transformation of the nucleotide frequency params in `phi`;
            all entries are > 0 and < 1.
        `freeparams` (list of strings)
            List of the model parameters that are free to be optimized. Can
            include `eta` but **not** `phi`, as they represent same thing.
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
        `ln_pi_codon` (`numpy.ndarray` floats, shape `(nsites, N_CODON)`)
            Natural logarithm of `pi_codon`.
        `piAx_piAy` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
            `piAx_piAy[r][x][y]` is ratio of preference for amino acid
            encoded by codon `x` divided by pref for that encoded by `y`.
        `piAx_piAy_beta` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
            Equal to `piAx_piAy` raised to the power of `beta`.
        `dPrxy` (dict)
            Keyed by each string in `freeparams`, each value is `numpy.ndarray`
            of floats giving derivative of `Prxy` with respect to that parameter.
            The shape of each array `(nsites, N_CODON, N_CODON)` except for 
            *eta*, for which it is `(N_NT - 1, nsites, N_CODON, N_CODON)` with
            the first index ranging over each element in `eta`.
        `dprx` (dict)
            Keyed by each string in `freeparams`, each value is `numpy.ndarray`
            of floats giving derivative of `prx` with respect to that parameter,
            or 0 if if `prx` does not depend on parameter. The shape of each
            array is `(nsites, N_CODON)` except for *eta*, for which it is
            `(N_NT - 1, nsites, N_CODON)` with the first index ranging over
            each element in `eta`.
    """

    def __init__(self, prefs, kappa=1.0, omega=1.0, beta=1.0, phi=scipy.ones(N_NT) / N_NT,
            freeparams=['kappa', 'omega', 'beta', 'eta']):
        """Initialize an `ExpCM` object.

        Args: 
            `prefs` (list)
                List of dicts giving amino-acid preferences for
                each site. Each dict keyed by amino acid letter
                codes, value is pref > 0 and < 1. Must sum to 1 
                at each site. All values are adjusted to be >=
                ALMOST_ZERO.
            `kappa`, `omega`, `beta`, `phi`, `freeparams`.
                Parameters with meanings described in main class doc string.
        """
        self.nsites = len(prefs)
        assert self.nsites > 0, "No preferences specified"

        allowedfreeparams = ['kappa', 'omega', 'beta', 'eta']
        assert all(map(lambda x: x in allowedfreeparams, freeparams)),\
                "Invalid entry in freeparams\nGot: {0}\nAllowed: {1}".format(
                ', '.join(freeparams), ', '.join(allowedparams))
        self.freeparams = list(freeparams)

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
        self.ln_pi_codon = scipy.full((self.nsites, N_CODON), -1, dtype='float')
        self._update_ln_pi_codon()
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
        self.dPrxy = {}
        self.dprx = {}
        for param in self.freeparams:
            if param in ['kappa', 'omega', 'beta']:
                self.dPrxy[param] = scipy.zeros((self.nsites, N_CODON, N_CODON), 
                        dtype='float')
                if param == 'beta':
                    self.dprx[param] = scipy.zeros((self.nsites, N_CODON), dtype='float')
                else:
                    self.dprx[param] = 0.0
            elif param == 'eta':
                self.dPrxy[param] = scipy.zeros((N_NT - 1, self.nsites, N_CODON,
                        N_CODON), dtype='float')
                self.dprx[param] = scipy.zeros((N_NT - 1, self.nsites, N_CODON),
                        dtype='float')
            else:
                raise ValueError("Unrecognized param {0}".format(param))

        # indexes diagonals in square matrices
        self._diag_indices = scipy.diag_indices(N_CODON)
        
        self.updateParams({}, update_all=True)

    def updateParams(self, newvalues, update_all=False):
        """Update model params.

        This method is the **only** acceptable way to update `ExpCM`
        model parameters. This method automatically updates
        any other dependent `ExpCM` attributes after updating
        the parameters to the new value.

        Args:
            `newvalues` (dict)
                Can be keyed by any parameter in `freeparams`.
            `update_all` (bool)
                If `True`, update all dependent attributes using
                current values of model parameters.

        """
        assert all(map(lambda x: x in self.freeparams, newvalues.keys())),\
                "Invalid entry in newvalues: {0}\nfreeparams: {1}".format(
                ', '.join(newvalues.keys()), ', '.join(self.freeparams))
        changed = set([]) # contains string names of changed params

        if ('eta' in newvalues) and (newvalues['eta'] != self.eta).any():
            self.eta = newvalues['eta']
            assert (isinstance(self.eta, scipy.ndarray) and self.eta.shape == (3,) 
                    and self.eta.dtype == float), "eta not array of 3 floats"
            assert (self.eta > 0).all() and (self.eta < 1).all(),\
                    "eta must be > 0 and < 1"
            changed.add('eta')

        if ('kappa' in newvalues) and newvalues['kappa'] != self.kappa:
            self.kappa = newvalues['kappa']
            assert isinstance(self.kappa, float) and self.kappa > 0,\
                    "kappa must be float > 0"
            changed.add('kappa')

        if ('omega' in newvalues) and newvalues['omega'] != self.omega:
            self.omega = newvalues['omega']
            assert isinstance(self.omega, float) and self.omega > 0,\
                    "omega must be float > 0"
            changed.add('omega')

        if ('beta' in newvalues) and newvalues['beta'] != self.beta:
            self.beta = newvalues['beta']
            assert isinstance(self.beta, float) and self.beta > 0,\
                    "beta must be float > 0"
            changed.add('beta')

        if update_all or ('eta' in changed):
            self._update_phi()
            self._update_Qxy()
            self._update_qx()
        elif 'kappa' in changed:
            self._update_Qxy()

        if update_all or ('beta' in changed):
            self._update_piAx_piAy_beta()
            self._update_Frxy()
            self._update_frx()
        elif 'omega' in changed:
            self._update_Frxy()

        if update_all or ('beta' in changed) or ('eta' in changed):
            self._update_prx()
            self._update_dprx()

        if update_all or changed:
            self._update_Prxy()
            self._update_dPrxy()

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

    def _update_ln_pi_codon(self):
        """Update `ln_pi_codon` using current `pi_codon`."""
        with scipy.errstate(divide='raise', under='raise', over='raise',
                invalid='raise'):
            self.ln_pi_codon = scipy.log(self.pi_codon)

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
        scipy.copyto(self.Frxy, self.omega, where=(scipy.logical_and(CODON_NONSYN,
                1 == self.piAx_piAy_beta)))
        assert not scipy.isnan(self.Frxy).any(), "NaN values in Frxy"
        assert not scipy.isinf(self.Frxy).any(), "Infinite values in Frxy"

    def _update_frx(self):
        """Update `frx` using current `pi_codon` and `beta`."""
        scipy.copyto(self.frx, self.pi_codon**self.beta)

    def _update_Prxy(self):
        """Update `Prxy` using current `Frxy` and `Qxy`."""
        scipy.copyto(self.Prxy, self.Frxy * self.Qxy)
        self._fill_diagonals(self.Prxy)

    def _update_prx(self):
        """Update `prx` using current `frx` and `qx`."""
        self.prx = self.frx * self.qx
        with scipy.errstate(divide='raise', under='raise', over='raise',
                invalid='raise'):
            for r in range(self.nsites):
                self.prx[r] /= self.prx[r].sum()

    def _update_dPrxy(self):
        """Update `dPrxy`."""
        if 'kappa' in self.freeparams:
            scipy.copyto(self.dPrxy['kappa'], self.Prxy / self.kappa, 
                    where=CODON_TRANSITION)
            self._fill_diagonals(self.dPrxy['kappa'])
        if 'omega' in self.freeparams:
            scipy.copyto(self.dPrxy['omega'], self.Prxy / self.omega, 
                    where=CODON_NONSYN)
            self._fill_diagonals(self.dPrxy['omega'])
        if 'beta' in self.freeparams:
            scipy.copyto(self.dPrxy['beta'], (self.Prxy * (1 - self.Frxy / self.omega
                    * self.piAx_piAy_beta)) / self.beta, where=CODON_NONSYN)
            self._fill_diagonals(self.dPrxy['beta'])
        if 'eta' in self.freeparams:
            for i in range(N_NT - 1):
                for w in range(i, N_NT):
                    scipy.copyto(self.dPrxy['eta'][i], self.Prxy / (self.eta[i]
                            - int(i == w)), where=CODON_NT_MUT[w])
                self._fill_diagonals(self.dPrxy['eta'][i])
            
    def _update_dprx(self):
        """Update `dprx`."""
        if 'beta' in self.freeparams:
            for r in range(self.nsites):
                scipy.copyto(self.dprx['beta'][r], self.prx[r] * (self.ln_pi_codon[r]
                        - scipy.dot(self.ln_pi_codon[r], self.prx[r])))
        if 'eta' in self.freeparams:
            boolterm = scipy.ndarray(N_CODON, dtype='float')
            with scipy.errstate(divide='raise', under='raise', over='raise',
                    invalid='raise'):
                for i in range(N_NT - 1):
                    boolterm.fill(0)
                    for j in range(3):
                        boolterm += ((i <= CODON_NT_INDEX[j]).astype('float') / 
                                (self.eta[i] - (i == CODON_NT_INDEX[j]).astype(
                                'float')))
                    for r in range(self.nsites):
                        scipy.copyto(self.dprx['eta'][i][r], self.prx[r] * (boolterm
                                - scipy.dot(boolterm, self.prx[r]) / self.prx[r].sum()))

    def _fill_diagonals(self, m):
        """Fills diagonals of `nsites` matrices in `m` so rows sum to 0."""
        assert m.shape == (self.nsites, N_CODON, N_CODON)
        for r in range(self.nsites):
            scipy.fill_diagonal(m[r], 0)
            m[r][self._diag_indices] = -scipy.sum(m[r], axis=1)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
