"""Substitution models.

All models defined here are subclasses of abstract base class `Model`.

For all models, nucleotides, amino acids, and codons are indexed
by integers 0, 1, ... using the indexing schemes defined
in `phydmslib.constants`.
"""


import math
import six
import abc
import scipy
import scipy.linalg
from phydmslib.constants import *


class Model(six.with_metaclass(abc.ABCMeta)):
    """Substitution model abstract base class.
    
    Specifies required methods / attributes of substitution models.
    """

    @abc.abstractmethod
    def M(self, t):
        """Matrix exponential `M(mu * t) = exp(mu * t * P)`.
       
        Args:
            `t` (float > 0)
                The branch length.

        Returns:
            `M` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
                `M(t)[r][x][y]` is probability `r` changes from `x` to `y`
                in time `t`.
        """
        pass

    @abc.abstractmethod
    def dM(self, t, param):
        """Derivative of `M(t)` with respect to `param`.

        Args:
            `t` (float > 0)
                The branch length.
            `param` (string in `freeparams`)
                Differentiate with respect to this model parameter.

        Returns:
            `dM_param` (`numpy.ndarray` floats)
                If `param` is a float, then `dM_param[r][x][y]` is 
                derivative of `M(t)[r][x][y]` with respect to `param`.
                If `param` is an array, then `dM_param[i][r][x][y]`
                is derivative of `M(t)[r][x][y]` with respect to
                `param[i]`.
        """
        pass

    @abc.abstractmethod
    def updateParams(self, newvalues, update_all=False):
        """Update model params.

        This method is the **only** acceptable way to update 
        model parameters; other dependent model parameters are
        automatically updated as needed if you use this method.

        Args:
            `newvalues` (dict)
                Can be keyed by any parameter in `freeparams`.
            `update_all` (bool)
                If `True`, update all dependent attributes using
                current values of model parameters.
        """
        pass

    @abc.abstractproperty
    def nsites(self):
        """An `int` giving the number of codon sites."""
        pass

    @abc.abstractproperty
    def freeparams(self):
        """Model parameters to be optimized (a list of strings)."""
        pass

    @abc.abstractproperty
    def mu(self):
        """Mutation rate, scales branch lengths in `M` and `dM`."""
        pass

    @abc.abstractproperty
    def ALLOWEDPARAMS(self):
        """List of all strings that can be included in `freeparam`."""
        pass

    @abc.abstractproperty
    def PARAMLIMITS(self):
        """Dict of tuples giving lower and upper allowed param values.

        For each `param` in `ALLOWEDPARAMS`, `PARAMLIMITS[param]`
        is a 2-tuple of floats giving upper and lower allowed
        values of `param` (or of each element in `param` if
        `param` is an array).
        """
        pass



class ExpCM(Model):
    """Experimentally informed codon models (ExpCM) for a gene.

    See `__init__` method for how to initialize an `ExpCM`.

    Attributes should **only** be updated via the `updateParams`
    method, do **not** set attributes directly.

    Attributes of `ExpCM` instances: 
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
            computed from `eta`. We use `eta` rather than `phi`
            as the free parameter during optimization.
        `eta` (`numpy.ndarray` of floats, length `N_NT - 1`)
            Transformation of the nucleotide frequency params in `phi`;
            all entries are > 0 and < 1. We use `eta` rather than `phi`
            as the free parameter during optimization.
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
            Keyed by strings in `freeparams`, each value is `numpy.ndarray`
            of floats giving derivative of `prx` with respect to that param,
            or 0 if if `prx` does not depend on parameter. The shape of each
            array is `(nsites, N_CODON)` except for *eta*, for which it is
            `(N_NT - 1, nsites, N_CODON)` with the first index ranging over
            each element in `eta`.
        `D` (`numpy.ndarray` floats, shape `(nsites, N_CODON)`)
            Element `r` is diagonal of :math:`\mathbf{D_r}` diagonal matrix where
            :math:`\mathbf{P_r} = \mathbf{A_r}^{-1} \mathbf{D_r} \mathbf{A_r}`
        `Dxx` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
            Element `[r][x][y]` is `D[r][x]`.
        `Dyy` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
            Element `[r][x][y]` is `D[r][y]`.
        `Dxx_Dyy` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
            Equals `Dxx - Dyy`.
        `A` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
            Element r is :math:`\mathbf{A_r}` matrix.
        `Ainv` (`numpy.ndarray` floats, shape `(nsites, N_CODON, N_CODON)`)
            Element r is :math:`\mathbf{A_r}^{-1}` matrix.
        `B` (dict)
            Keyed by strings in `freeparams`, each value is `numpy.ndarray`
            of floats giving `B` matrix for that parameter. The shape
            is `(nsites, N_CODON, N_CODON)` except for *eta*, for which 
            it is `(N_NT - 1, nsites, N_CODON, N_CODON)` with the first
            index ranging over each element in `eta`.
    """

    # class variables
    _ALLOWEDPARAMS = ['kappa', 'omega', 'beta', 'eta', 'mu']
    _PARAMLIMITS = {'kappa':(0.01, 100.0),
                   'omega':(0.01, 100.0),
                   'beta':(0.01, 5.0),
                   'eta':(0.01, 0.99),
                   'phi':(0.001, 0.999),
                   'pi':(0.002, 0.998),
                   'mu':(1.0e-3, 1.0e3),
                  }

    @property
    def ALLOWEDPARAMS(self):
        """See docs for `Model` abstract base class."""
        return self._ALLOWEDPARAMS

    @property
    def PARAMLIMITS(self):
        """See docs for `Model` abstract base class."""
        return self._PARAMLIMITS

    def checkParam(self, param, value):
        """Checks if `value` is allowable value for `param`.

        Raises except if `value` is not acceptable, otherwise
        returns `None` if value is acceptable.
        """
        assert param in self.PARAMLIMITS, "Invalid param: {0}".format(param)
        (lowlim, highlim) = self.PARAMLIMITS[param]
        if param in ['kappa', 'omega', 'beta', 'pi', 'mu']:
            if not isinstance(value, float):
                raise ValueError("{0} must be a float, not a {1}".format(
                        param, type(value)))
            if not (lowlim <= value <= highlim):
                raise ValueError("{0} must be >= {1} and <= {2}, not {3}".format(
                        param, lowlim, highlim, value))
        elif param in ['eta', 'phi']:
            shape = {'eta':(N_NT - 1,),
                     'phi':(N_NT,),
                    }
            if not (isinstance(value, scipy.ndarray) and value.shape == shape[param]
                    and value.dtype == 'float'):
                raise ValueError("{0} must be ndarray floats, shape {1}, " + 
                        "not {2}".format(param, shape, value))
            if not ((lowlim <= value).all() and (value <= highlim).all()):
                raise ValueError("{0} must be >= {1} and <= {2}, not {3}".format(
                        param, lowlim, highlim, value))
        else:
            raise ValueError("Can't handle {0}".format(param))

    def __init__(self, prefs, kappa=1.0, omega=1.0, beta=1.0, mu=1.0,
            phi=scipy.ones(N_NT) / N_NT,
            freeparams=['kappa', 'omega', 'beta', 'mu', 'eta']):
        """Initialize an `ExpCM` object.

        Args: 
            `prefs` (list)
                List of dicts giving amino-acid preferences for
                each site. Each dict keyed by amino acid letter
                codes, value is pref > 0 and < 1. Must sum to 1 
                at each site. 
            `kappa`, `omega`, `beta`, `mu`, `phi`
                Model params described in main class doc string.
            `freeparams` (list of strings)
                Specifies free parameters.
        """
        self._nsites = len(prefs)
        assert self.nsites > 0, "No preferences specified"

        assert all(map(lambda x: x in self.ALLOWEDPARAMS, freeparams)),\
                "Invalid entry in freeparams\nGot: {0}\nAllowed: {1}".format(
                ', '.join(freeparams), ', '.join(self.ALLOWEDPARAMS))
        self._freeparams = list(freeparams) # underscore as `freeparams` is property

        # put prefs in pi
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
                self.checkParam('pi', prefs[r][aa])
                self.pi[r][a] = prefs[r][aa]
            self.pi[r] /= self.pi[r].sum() # renormalize to sum to one 

        # set up attributes defined solely in terms of preferences
        self.pi_codon = scipy.full((self.nsites, N_CODON), -1, dtype='float')
        self.ln_pi_codon = scipy.full((self.nsites, N_CODON), -1, dtype='float')
        self.piAx_piAy = scipy.full((self.nsites, N_CODON, N_CODON), -1, 
                dtype='float')
        self._update_pi_vars()

        # construct eta from phi 
        self.checkParam('phi', phi)
        assert abs(1 - phi.sum()) <= ALMOST_ZERO, "phi doesn't sum to 1"
        self.phi = phi.copy()
        self.phi /= self.phi.sum()
        eta = scipy.ndarray(N_NT - 1, dtype='float')
        etaprod = 1.0
        for w in range(N_NT - 1):
            eta[w] = 1.0 - self.phi[w] / etaprod
            etaprod *= eta[w]

        # set attributes to calling params
        self._mu = mu # underscore as `mu` is property
        self.kappa = kappa
        self.omega = omega
        self.beta = beta
        self.eta = eta
        for (name, value) in [('kappa', self.kappa), ('omega', self.omega),
                ('beta', self.beta), ('eta', self.eta), ('mu', self.mu)]:
            self.checkParam(name, value)

        # define other params, initialized appropriately
        self.piAx_piAy_beta = scipy.zeros((self.nsites, N_CODON, N_CODON), 
                dtype='float')
        self.Prxy = scipy.zeros((self.nsites, N_CODON, N_CODON), dtype='float')
        self.prx = scipy.zeros((self.nsites, N_CODON), dtype='float')
        self.Qxy = scipy.zeros((N_CODON, N_CODON), dtype='float')
        self.qx = scipy.zeros(N_CODON, dtype='float')
        self.Frxy = scipy.ones((self.nsites, N_CODON, N_CODON), dtype='float')
        self.frx = scipy.zeros((self.nsites, N_CODON), dtype='float')
        self.D = scipy.zeros((self.nsites, N_CODON), dtype='float')
        self.Dxx = scipy.zeros((self.nsites, N_CODON, N_CODON), dtype='float')
        self.Dyy = scipy.zeros((self.nsites, N_CODON, N_CODON), dtype='float')
        self.Dxx_Dyy = scipy.zeros((self.nsites, N_CODON, N_CODON), dtype='float')
        self.A = scipy.zeros((self.nsites, N_CODON, N_CODON), dtype='float')
        self.Ainv = scipy.zeros((self.nsites, N_CODON, N_CODON), dtype='float')
        self.dPrxy = {}
        self.B = {}
        self.dprx = {}
        for param in self.freeparams:
            if param in ['kappa', 'omega', 'beta']:
                self.dPrxy[param] = scipy.zeros((self.nsites, N_CODON, N_CODON), 
                        dtype='float')
                self.B[param] = scipy.zeros((self.nsites, N_CODON, N_CODON), 
                        dtype='float')
                if param == 'beta':
                    self.dprx[param] = scipy.zeros((self.nsites, N_CODON), dtype='float')
                else:
                    self.dprx[param] = 0.0
            elif param == 'eta':
                self.dPrxy[param] = scipy.zeros((N_NT - 1, self.nsites, N_CODON,
                        N_CODON), dtype='float')
                self.B[param] = scipy.zeros((N_NT - 1, self.nsites, N_CODON,
                        N_CODON), dtype='float')
                self.dprx[param] = scipy.zeros((N_NT - 1, self.nsites, N_CODON),
                        dtype='float')
            elif param == 'mu':
                pass 
            else:
                raise ValueError("Unrecognized param {0}".format(param))

        # indexes diagonals in square matrices
        self._diag_indices = scipy.diag_indices(N_CODON)

        self._cached_M = {} # caches results of calls to M
        self._cached_dM = {} # caches results of calls to dM
        self.updateParams({}, update_all=True)

    @property
    def nsites(self):
        """See docs for `Model` abstract base class."""
        return self._nsites

    @property
    def freeparams(self):
        """See docs for `Model` abstract base class."""
        return self._freeparams

    @property
    def mu(self):
        """See docs for `Model` abstract base class."""
        return self._mu

    def updateParams(self, newvalues, update_all=False):
        """See docs for `Model` abstract base class."""
        assert all(map(lambda x: x in self.freeparams, newvalues.keys())),\
                "Invalid entry in newvalues: {0}\nfreeparams: {1}".format(
                ', '.join(newvalues.keys()), ', '.join(self.freeparams))
        changed = set([]) # contains string names of changed params

        for (name, value) in newvalues.items():
            self.checkParam(name, value)

        if ('mu' in newvalues) and (newvalues['mu'] != self.mu):
            self._mu = newvalues['mu']
            changed.add('mu')

        if ('eta' in newvalues) and (newvalues['eta'] != self.eta).any():
            self.eta = newvalues['eta']
            changed.add('eta')

        if ('kappa' in newvalues) and newvalues['kappa'] != self.kappa:
            self.kappa = newvalues['kappa']
            changed.add('kappa')

        if ('omega' in newvalues) and newvalues['omega'] != self.omega:
            self.omega = newvalues['omega']
            changed.add('omega')

        if ('beta' in newvalues) and newvalues['beta'] != self.beta:
            self.beta = newvalues['beta']
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

        if update_all or (changed and changed != set(['mu'])):
            self._cached_M = {}
            self._cached_dM = {}
            self._update_Prxy()
            self._update_Prxy_diag()
            self._update_dPrxy()
            self._update_B()
        elif changed: # only changed mu
            self._cached_M = {}
            self._cached_dM = {}

    def M(self, t):
        """See docs for method in `Model` abstract base class."""
        if t in self._cached_M:
            return self._cached_M[t]
        assert isinstance(t, float) and t > 0, "Invalid t: {0}".format(t)
        # swap axes commands allow broadcasting to multiply D like diagonal matrix
        M = scipy.matmul((self.A.swapaxes(0, 1) * scipy.exp(self.D * self.mu * t)
                ).swapaxes(1, 0), self.Ainv)
        self._cached_M[t] = M
        return M

    def dM(self, t, param):
        """See docs for method in `Model` abstract base class."""
        key = (t, param)
        if key in self._cached_dM:
            return self._cached_dM[key]
        assert isinstance(t, float) and t > 0, "Invalid t: {0}".format(t)
        assert param in self.freeparams, "Invalid param: {0}".format(param)
        if param == 'mu':
            dM_param = t * scipy.matmul(self.Prxy, self.M(t))
            self._cached_dM[key] = dM_param
            return dM_param
        mut = self.mu * t
        with scipy.errstate(divide='raise', under='ignore', over='raise',
                invalid='ignore'):
            V = (scipy.exp(mut * self.Dxx) - scipy.exp(mut * self.Dyy)) / self.Dxx_Dyy
        scipy.copyto(V, mut * scipy.exp(mut * self.Dxx), where=
                scipy.fabs(self.Dxx_Dyy) < ALMOST_ZERO)
        dM_param = scipy.matmul(self.A, scipy.matmul(self.B[param] * V, self.Ainv))
        self._cached_dM[key] = dM_param
        return dM_param

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

    def _update_pi_vars(self):
        """Update `pi_codon`, `ln_pi_codon`, `piAx_piAy` from current `pi`."""
        with scipy.errstate(divide='raise', under='raise', over='raise',
                invalid='raise'):
            for r in range(self.nsites):
                for a in range(N_AA):
                    scipy.copyto(self.pi_codon[r], self.pi[r][a], 
                    where=(CODON_TO_AA == a))
                pim = scipy.tile(self.pi_codon[r], (N_CODON, 1)) # [x][y] is piAy
                scipy.copyto(self.piAx_piAy[r], pim.transpose() / pim)
            self.ln_pi_codon = scipy.log(self.pi_codon)

    def _update_piAx_piAy_beta(self):
        """Update `piAx_piAy_beta` from `piAx_piAy` and `beta`."""
        with scipy.errstate(divide='raise', under='raise', over='raise', 
                invalid='raise'):
            scipy.copyto(self.piAx_piAy_beta, self.piAx_piAy**self.beta)

    def _update_Frxy(self):
        """Update `Frxy` from `piAx_piAy_beta`, `omega`, and `beta`."""
        with scipy.errstate(divide='raise', under='raise', over='raise', 
                        invalid='ignore'):
            scipy.copyto(self.Frxy, -self.omega * scipy.log(self.piAx_piAy_beta)
                    / (1 - self.piAx_piAy_beta), where=CODON_NONSYN)
        scipy.copyto(self.Frxy, self.omega, where=(scipy.logical_and(CODON_NONSYN,
                scipy.fabs(1 - self.piAx_piAy_beta) < ALMOST_ZERO)))

    def _update_frx(self):
        """Update `frx` using current `pi_codon` and `beta`."""
        scipy.copyto(self.frx, self.pi_codon**self.beta)

    def _update_Prxy(self):
        """Update `Prxy` using current `Frxy` and `Qxy`."""
        scipy.copyto(self.Prxy, self.Frxy * self.Qxy)
        self._fill_diagonals(self.Prxy)

    def _update_Prxy_diag(self):
        """Updates `D`, `A`, `Ainv`, `Dxx`, `Dyy`, `Dxx_Dyy` from `Prxy`, `prx`."""
        for r in range(self.nsites):
            pr_half = self.prx[r]**0.5
            pr_neghalf = self.prx[r]**-0.5
            #symm_pr = scipy.dot(scipy.diag(pr_half), scipy.dot(self.Prxy[r], scipy.diag(pr_neghalf)))
            symm_pr = (pr_half * (self.Prxy[r] * pr_neghalf).transpose()).transpose()
            #assert scipy.allclose(symm_pr, symm_pr.transpose())
            (evals, evecs) = scipy.linalg.eigh(symm_pr)
            #assert scipy.allclose(scipy.linalg.inv(evecs), evecs.transpose())
            #assert scipy.allclose(symm_pr, scipy.dot(evecs, scipy.dot(scipy.diag(evals), evecs.transpose())))
            scipy.copyto(self.D[r], evals)
            scipy.copyto(self.Ainv[r], evecs.transpose() * pr_half)
            scipy.copyto(self.A[r], (pr_neghalf * evecs.transpose()).transpose())
            scipy.copyto(self.Dxx[r], scipy.tile(self.D[r], (N_CODON, 1)).transpose())
            scipy.copyto(self.Dyy[r], scipy.tile(self.D[r], (N_CODON, 1)))
        scipy.copyto(self.Dxx_Dyy, self.Dxx - self.Dyy)

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

    def _update_B(self):
        """Update `B`."""
        for param in self.freeparams:
            if param != 'mu':
                scipy.copyto(self.B[param], scipy.matmul(self.Ainv, scipy.matmul(
                        self.dPrxy[param], self.A)))

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
