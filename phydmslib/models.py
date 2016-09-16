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

    Attributes: do **not** modify except via `updateParams` method.
        `nsites` (int)
            Number of codon sites in gene.
        `pi` (`numpy.ndarray` of floats, shape `(nsites, N_AA)`)
            `pi[r][x]` gives amino-acid preference of site `r`
            for amino-acid `x`.
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
    """

    def __init__(self, prefs, kappa=1.0, omega=1.0, beta=1.0,
            phi=scipy.ones(N_NT) / N_NT):
        """Initialize an `ExpCM` object.

        Args: 
            `prefs` (list)
                List of dicts giving amino-acid preferences for
                each site. Each dict keyed by amino acid letter
                codes, value is pref. Must sum to 1 at each site.
            `kappa`, `omega`, `beta`, and `phi`
                Meanings described in main class doc string.
        """
        self.nsites = len(prefs)
        assert self.nsites > 0, "No preferences specified"

        # put prefs in pi, adjusting to be > ALMOST_ZERO
        self.pi = scipy.ndarray((self.nsites, N_AA))
        assert (isinstance(prefs, list) and 
                all([isinstance(x, dict) for x in prefs])),\
                "prefs is not a list of dicts"
        for r in range(self.nsites):
            assert set(prefs[r].keys()) == set(AA_TO_INDEX.keys()),\
                    "prefs not keyed by amino acids for site {0}".format(r)
            assert abs(1 - sum(prefs[r].values())) <= ALMOST_ZERO,\
                    "prefs don't sum to one for site {0}".format(r)
            assert all([0 <= pix <= 1 for pix in prefs[r].values()]),\
                    "prefs aren't all >= 0 and <= 1 for site {0}".format(r)
            for (x, aa) in INDEX_TO_AA.items():
                self.pi[r][x] = max(ALMOST_ZERO, prefs[r][aa])
            self.pi[r] /= self.pi[r].sum() # renormalize to sum to one

        # construct eta from phi after ensuring > ALMOST_ZERO
        assert len(phi) == N_NT, "not the right number of entries in phi"
        assert abs(1 - phi.sum()) <= ALMOST_ZERO, "phi doesn't sum to 1"
        assert (phi > 0).all() and (phi < 1).all(),\
                "phi entries not > 0 and < 1"
        self.phi = phi.copy()
        self.phi[self.phi < ALMOST_ZERO] = ALMOST_ZERO
        self.phi /= self.phi.sum()
        eta = scipy.ndarray(N_NT - 1)
        etaprod = 1.0
        for w in range(N_NT - 1):
            eta[w] = 1.0 - self.phi[w] / etaprod
            etaprod *= eta[w]

        # now set the calling params as object attributes
        assert kappa != None
        assert omega != None
        assert beta != None
        self.kappa = self.omega = self.beta = self.eta = None
        self.updateParams(kappa=kappa, omega=omega, beta=beta, eta=eta)


    def updateParams(self, kappa=None, omega=None, beta=None, eta=None):
        """Update model params.

        This method is the **only** acceptable way to update `ExpCM`
        model parameters. Any arguments that are not set to the
        default value of `None` are updated if they changed. This method
        then automatically updates any other `ExpCM` attributes
        that need to be changed as a result of the parameter update.

        Args:
            `kappa`, `omega`, `beta`, and `phi`
                Meanings described in main class doc string.
        """
        if eta is not None and eta != self.eta:
            assert (isinstance(eta, scipy.ndarray) and eta.shape == (3,) 
                    and eta.dtype == float), "eta not array of 3 floats"
            assert (eta > 0).all() and (eta < 1).all(),\
                    "eta must be > 0 and < 1"
            self.eta = eta
            # update self.phi
            etaprod = 1.0
            for w in range(N_NT - 1):
                self.phi[w] = etaprod * (1 - self.eta[w])
                etaprod *= self.eta[w]
            self.phi[N_NT - 1] = etaprod

        if kappa != None and kappa != self.kappa:
            assert isinstance(kappa, float) and kappa > 0,\
                    "kappa must be float > 0"
            self.kappa = kappa

        if omega != None and omega != self.omega:
            assert isinstance(omega, float) and omega > 0,\
                    "omega must be float > 0"
            self.omega = omega

        if beta != None and beta != self.beta:
            assert isinstance(beta, float) and beta > 0,\
                    "beta must be float > 0"
            self.beta = beta



if __name__ == '__main__':
    import doctest
    doctest.testmod()
