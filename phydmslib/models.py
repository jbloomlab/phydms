"""Substitution models.

All models defined here are subclasses of abstract base class `Model`.

Nucleotides, amino acids, and codons are indexed by integers 0, 1, ...
using the indexing schemes defined in `phydmslib.constants`.
"""


import math
import copy
import functools
import six
import abc
import scipy
import scipy.misc
import scipy.optimize
import scipy.linalg
import scipy.stats
import scipy.special
from phydmslib.numutils import *
from phydmslib.constants import *


class Model(six.with_metaclass(abc.ABCMeta)):
    """Substitution model abstract base class.

    Specifies required methods / attributes of substitution models.
    """

    @abc.abstractproperty
    def stationarystate(self):
        """Stationary state of substitution model.

        A `numpy.ndarray` of floats, shape `(nsites, N_CODON)`.
        Element `stationarystate[r][x]` is the stationary
        state probability of codon `x` at site `r`.
        """
        pass

    @abc.abstractmethod
    def dstationarystate(self, param):
        """Derivative of `stationarystate` with respect to `param`.

            Args:
                `param` (string in `freeparams`)

            Returns:
                `dstationarystate` (`numpy.ndarray` of floats)
                    If `param` is a float, then `dstationarystate[r][x]`
                    is derivative of `stationarystate[r][x]` with respect
                    to `param`. If `param` is an array, then
                    `dstationarystate[i][r][x]` is derivative of
                    `stationarystate[r][x]` with respect to `param[i]`.
        """
        pass

    @abc.abstractproperty
    def paramsReport(self):
        """Reports current values of independent model parameters.

        Returns a dictionary keyed by parameter name with value being the 
        parameter value. Does **not** include `mu` as this is confounded
        with branch lengths.
        """
        pass

    @abc.abstractproperty
    def branchScale(self):
        """Factor to scale branch lengths to substitutions per site.

        When we use a branch length of `t = 1`, how many substitutions
        per site do we expect averaged over all sites? This method
        returns that number. So to scale branch lengths to units of
        substitution per site, multiply by this number.
        """
        pass

    @abc.abstractmethod
    def M(self, t, tips=None, gaps=None):
        """Matrix exponential `M(mu * t) = exp(mu * t * P)`.

        If we are considering branch to tip node, you can use
        `tips` and `gaps` to get only the column needed for
        likelihood calculation. This is computationally cheaper.

        Args:
            `t` (float > 0)
                The branch length.
            `tips` (`None` or `numpy.ndarray` of `int`, length `nsites`)
                If `None`, return the full matrix exponential.
                Otherwise, should be an array with element `r`
                giving the codon at tip node `r` (just put 0
                if codon is a gap and see `gaps` below).
            `gaps` (`None` or `numpy.ndarray` of `int`)
                Only meaningful if using `tips`. In this case,
                array should list `r` for all sites where tip
                node codon is a gap. `None` is equivalent to
                no gaps (also can be empty array).

        Returns:
            `M` (`numpy.ndarray` of floats)
                If not using `tips`, shape is `(nsites, N_CODON, N_CODON)`
                    with `M(t)[r][x][y]` probability `r` changes from `x`
                    to `y` in time `t`.
                If using `tips`, shape is `(nsites, N_CODON)` with
                    `M(t)[r][x] probability r changes from `x` to `tips[r]`
                    in time `t`. If `r` is in `gap`, then `M(t)[r][x]` is
                    one.
        """
        pass

    @abc.abstractmethod
    def dM(self, t, param, Mt, tips=None, gaps=None):
        """Derivative of `M(t)` with respect to `param`.

        Args:
            `t` (float > 0)
                The branch length.
            `param` (string in `freeparams`)
                Differentiate with respect to this model parameter.
            `Mt` (`numpy.ndarray`.)
                The current value of `M(t, tips, gaps)`. Typically this
                has already been pre-computed which is why it is passed
                here to avoid computing again. Otherwise pass
                `None` and it will be re-computed.
            `tips` and `gaps`
                Same meaning as for `M`.

        Returns:
            `dM_param` (`numpy.ndarray` floats)
                If `param` is a float, then `dM_param[r][x][y]`
                is derivative of `M(t)[r][x][y]` with respect to
                `param`. If `param` is an array, then
                `dM_param[i][r][x][y]` is derivative of `M(t)[r][x][y]`
                with respect to `param[i]`. If using `tips`, then
                `dM_param[r]` or `dM_params[i][r]` is only a
                single array providing the derivative of
                `M(t, tips, gaps)[r]` with respect to `param` or
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
        """List of all strings that can be included in `freeparams`."""
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

    Attributes (see also those inherited from `Model`):
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
            This attribute is equivalent to `stationarystate`.
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
            each element in `eta`. Equivalent to `dstationarystate[param]`.
        `D` (`numpy.ndarray` floats, shape `(nsites, N_CODON)`)
            Element `r` is diagonal of :math:`\mathbf{D_r}` diagonal matrix,
            :math:`\mathbf{P_r} = \mathbf{A_r}^{-1} \mathbf{D_r} \mathbf{A_r}`
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
    _REPORTPARAMS = ['kappa', 'omega', 'beta', 'phi']
    _PARAMLIMITS = {'kappa':(0.01, 100.0),
                   'omega':(0.01, 100.0),
                   'beta':(0.01, 10.0),
                   'eta':(0.01, 0.99),
                   'phi':(0.001, 0.999),
                   'pi':(0.002, 0.998),
                   'mu':(1.0e-3, 1.0e3),
                  }
    _PARAMTYPES = {'kappa':float,
                   'omega':float,
                   'beta':float,
                   'pi':float,
                   'mu':float,
                   'eta':(scipy.ndarray, (N_NT - 1,)),
                   'phi':(scipy.ndarray, (N_NT,)),
                  }

    @property
    def ALLOWEDPARAMS(self):
        """See docs for `Model` abstract base class."""
        return self._ALLOWEDPARAMS

    @property
    def PARAMLIMITS(self):
        """See docs for `Model` abstract base class."""
        return self._PARAMLIMITS

    def __init__(self, prefs, kappa=2.0, omega=0.5, beta=1.0, mu=1.0,
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
            for (a, aa) in INDEX_TO_AA.items():
                _checkParam('pi', prefs[r][aa], self.PARAMLIMITS, self._PARAMTYPES)
                self.pi[r][a] = prefs[r][aa]
            self.pi[r] /= self.pi[r].sum() # renormalize to sum to one

        # set up attributes defined solely in terms of preferences
        self.pi_codon = scipy.full((self.nsites, N_CODON), -1, dtype='float')
        self.ln_pi_codon = scipy.full((self.nsites, N_CODON), -1, dtype='float')
        self.piAx_piAy = scipy.full((self.nsites, N_CODON, N_CODON), -1,
                dtype='float')
        self._update_pi_vars()

        # construct eta from phi
        _checkParam('phi', phi, self.PARAMLIMITS, self._PARAMTYPES)
        assert abs(1 - phi.sum()) <= ALMOST_ZERO, "phi doesn't sum to 1"
        self.phi = phi.copy()
        self.phi /= self.phi.sum()
        self._eta_from_phi()

        # set attributes to calling params
        self._mu = mu # underscore as `mu` is property
        self.kappa = kappa
        self.omega = omega
        self.beta = beta
        for (name, value) in [('kappa', self.kappa), ('omega', self.omega),
                ('beta', self.beta), ('eta', self.eta), ('mu', self.mu)]:
            _checkParam(name, value, self.PARAMLIMITS, self._PARAMTYPES)

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
        self.A = scipy.zeros((self.nsites, N_CODON, N_CODON), dtype='float')
        self.Ainv = scipy.zeros((self.nsites, N_CODON, N_CODON), dtype='float')
        self.dPrxy = {}
        self.B = {}
        self.dprx = {}
        for param in self.freeparams:
            if param == 'eta':
                self.dPrxy[param] = scipy.zeros((N_NT - 1, self.nsites, N_CODON,
                        N_CODON), dtype='float')
                self.B[param] = scipy.zeros((N_NT - 1, self.nsites, N_CODON,
                        N_CODON), dtype='float')
                self.dprx[param] = scipy.zeros((N_NT - 1, self.nsites, N_CODON),
                        dtype='float')
            elif param == 'mu':
                self.dprx['mu'] = 0.0
            elif param in self._ALLOWEDPARAMS:
                self.dPrxy[param] = scipy.zeros((self.nsites, N_CODON, N_CODON),
                        dtype='float')
                self.B[param] = scipy.zeros((self.nsites, N_CODON, N_CODON),
                        dtype='float')
                self.dprx[param] = scipy.zeros((self.nsites, N_CODON), dtype='float')
            else:
                raise ValueError("Unrecognized param {0}".format(param))

        # indexes diagonals in square matrices
        self._diag_indices = scipy.diag_indices(N_CODON)

        self.updateParams({}, update_all=True)

    @property
    def stationarystate(self):
        """See docs for `Model` abstract base class."""
        return self.prx

    def dstationarystate(self, param):
        """See docs for `Model` abstract base class."""
        return self.dprx[param]

    @property
    def paramsReport(self):
        """See docs for `Model` abstract base class."""
        report = {}
        for param in self._REPORTPARAMS:
            pvalue = getattr(self, param)
            if isinstance(pvalue, float):
                report[param] = pvalue
            elif isinstance(pvalue, scipy.ndarray) and pvalue.shape == (N_NT,):
                for w in range(N_NT - 1):
                    report['{0}{1}'.format(param, INDEX_TO_NT[w])] = pvalue[w]
            else:
                raise RuntimeError("Unexpected param: {0}".format(param))
        return report

    @property
    def branchScale(self):
        """See docs for `Model` abstract base class."""
        bs = -(self.prx * scipy.diagonal(self.Prxy, axis1=1, axis2=2)
                ).sum() * self.mu / float(self.nsites)
        assert bs > 0
        return bs

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

    @mu.setter
    def mu(self, value):
        """Set new `mu` value."""
        self._mu = value

    def updateParams(self, newvalues, update_all=False):
        """See docs for `Model` abstract base class."""
        assert all(map(lambda x: x in self.freeparams, newvalues.keys())),\
                "Invalid entry in newvalues: {0}\nfreeparams: {1}".format(
                ', '.join(newvalues.keys()), ', '.join(self.freeparams))
        changed = set([]) # contains string names of changed params
        for (name, value) in newvalues.items():
            _checkParam(name, value, self.PARAMLIMITS, self._PARAMTYPES)
            if isinstance(value, scipy.ndarray):
                if (value != getattr(self, name)).any():
                    changed.add(name)
                    setattr(self, name, value)
            else:
                if value != getattr(self, name):
                    changed.add(name)
                    setattr(self, name, value)

        if update_all or changed:
            self._cached = {}

        # The order of the updating below is important.
        # If you change it, you may break either this class
        # **or** classes that inherit from it.
        # Note also that not all attributes need to be updated
        # for all possible parameter changes, but just doing it
        # this way is much simpler and adds negligible cost.
        if update_all or (changed and changed != set(['mu'])):
            self._update_piAx_piAy_beta()
            self._update_frx()
            self._update_phi()
            self._update_qx()
            self._update_prx()
            self._update_dprx()
            self._update_Qxy()
            self._update_Frxy()
            self._update_Prxy()
            self._update_Prxy_diag()
            self._update_dPrxy()
            self._update_B()

    def M(self, t, tips=None, gaps=None):
        """See docs for method in `Model` abstract base class."""
        assert isinstance(t, float) and t > 0, "Invalid t: {0}".format(t)
        with scipy.errstate(under='ignore'): # don't worry if some values 0
            if ('expD', t) not in self._cached:
                self._cached[('expD', t)] = scipy.exp(self.D * self.mu * t)
            expD = self._cached[('expD', t)]
            if tips is None:
                # swap axes to broadcast multiply D as diagonal matrix
                M = broadcastMatrixMultiply((self.A.swapaxes(0, 1) *
                        expD).swapaxes(1, 0), self.Ainv)
            else:
                M = broadcastMatrixVectorMultiply((self.A.swapaxes(0, 1)
                        * expD).swapaxes(1, 0), broadcastGetCols(
                        self.Ainv, tips))
                if gaps is not None:
                    M[gaps] = scipy.ones(N_CODON, dtype='float')
        return M

    def dM(self, t, param, Mt, tips=None, gaps=None):
        """See docs for method in `Model` abstract base class."""
        assert isinstance(t, float) and t > 0, "Invalid t: {0}".format(t)
        assert param in self.freeparams, "Invalid param: {0}".format(param)
        if param == 'mu':
            if tips is None:
                dM_param = broadcastMatrixMultiply(self.Prxy, Mt, alpha=t)
            else:
                dM_param = broadcastMatrixVectorMultiply(self.Prxy, Mt, alpha=t)
                if gaps is not None:
                    dM_param[gaps] = scipy.zeros(N_CODON, dtype='float')
            return dM_param

        paramval = getattr(self, param)
        if isinstance(paramval, float):
            paramisvec = False
        else:
            assert isinstance(paramval, numpy.ndarray) and paramval.ndim == 1
            paramisvec = True
            paramlength = paramval.shape[0]

        if ('expD', t) not in self._cached:
            self._cached[('expD', t)] = scipy.exp(self.D * self.mu * t)
        expD = self._cached[('expD', t)]

        if ('V', t) not in self._cached:
            if 'Dxx_Dyy' not in self._cached:
                Dyy = scipy.tile(self.D, (1, N_CODON)).reshape(
                        self.nsites, N_CODON, N_CODON)
                Dxx = scipy.array([Dyy[r].transpose() for r in
                        range(self.nsites)])
                self._cached['Dxx_Dyy'] = Dxx - Dyy
            Dxx_Dyy = self._cached['Dxx_Dyy']
            if 'Dxx_Dyy_lt_ALMOST_ZERO' not in self._cached:
                self._cached['Dxx_Dyy_lt_ALMOST_ZERO'] = scipy.fabs(
                        Dxx_Dyy) < ALMOST_ZERO
            Dxx_Dyy_lt_ALMOST_ZERO = self._cached['Dxx_Dyy_lt_ALMOST_ZERO']
            with scipy.errstate(divide='raise', under='ignore',
                    over='raise', invalid='ignore'):
                expDyy = scipy.tile(expD,(1, N_CODON)).reshape(
                        self.nsites, N_CODON, N_CODON)
                expDxx = scipy.array([expDyy[r].transpose() for r in
                        range(self.nsites)])
                V = (expDxx - expDyy) / Dxx_Dyy
            with scipy.errstate(under='ignore'): # OK if some values 0
                scipy.copyto(V, self.mu * t * expDxx, where=
                        Dxx_Dyy_lt_ALMOST_ZERO)
            self._cached[('V', t)] = V
        V = self._cached[('V', t)]

        with scipy.errstate(under='ignore'): # don't worry if some values 0
            if tips is None:
                if not paramisvec:
                    dM_param = broadcastMatrixMultiply(self.A,
                            broadcastMatrixMultiply(self.B[param]
                            * V, self.Ainv))
                else:
                    dM_param = scipy.ndarray((paramlength, self.nsites,
                            N_CODON, N_CODON), dtype='float')
                    for j in range(paramlength):
                        dM_param[j] = broadcastMatrixMultiply(self.A,
                                broadcastMatrixMultiply(self.B[param][j]
                                * V, self.Ainv))
            else:
                if not paramisvec:
                    dM_param = broadcastMatrixVectorMultiply(self.A,
                            broadcastGetCols(broadcastMatrixMultiply(
                            self.B[param] * V, self.Ainv), tips))
                else:
                    dM_param = scipy.ndarray((paramlength, self.nsites,
                            N_CODON), dtype='float')
                    for j in range(paramlength):
                        dM_param[j] = broadcastMatrixVectorMultiply(self.A,
                            broadcastGetCols(broadcastMatrixMultiply(
                            self.B[param][j] * V, self.Ainv), tips))
                if gaps is not None:
                    if not paramisvec:
                        dM_param[gaps] = scipy.zeros(N_CODON, dtype='float')
                    else:
                        dM_param[:, gaps] = scipy.zeros(N_CODON, dtype='float')
        return dM_param

    def _eta_from_phi(self):
        """Update `eta` using current `phi`."""
        self.eta = scipy.ndarray(N_NT - 1, dtype='float')
        etaprod = 1.0
        for w in range(N_NT - 1):
            self.eta[w] = 1.0 - self.phi[w] / etaprod
            etaprod *= self.eta[w]
        _checkParam('eta', self.eta, self.PARAMLIMITS, self._PARAMTYPES)

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
                self.piAx_piAy[r] = pim.transpose() / pim
            self.ln_pi_codon = scipy.log(self.pi_codon)

    def _update_piAx_piAy_beta(self):
        """Update `piAx_piAy_beta` from `piAx_piAy` and `beta`."""
        with scipy.errstate(divide='raise', under='raise', over='raise',
                invalid='raise'):
            self.piAx_piAy_beta = self.piAx_piAy**self.beta

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
        self.frx = self.pi_codon**self.beta

    def _update_Prxy(self):
        """Update `Prxy` using current `Frxy` and `Qxy`."""
        self.Prxy = self.Frxy * self.Qxy
        _fill_diagonals(self.Prxy, self._diag_indices)

    def _update_Prxy_diag(self):
        """Update `D`, `A`, `Ainv` from `Prxy`, `prx`."""
        for r in range(self.nsites):
            pr_half = self.prx[r]**0.5
            pr_neghalf = self.prx[r]**-0.5
            #symm_pr = scipy.dot(scipy.diag(pr_half), scipy.dot(self.Prxy[r], scipy.diag(pr_neghalf)))
            symm_pr = (pr_half * (self.Prxy[r] * pr_neghalf).transpose()).transpose()
            #assert scipy.allclose(symm_pr, symm_pr.transpose())
            (evals, evecs) = scipy.linalg.eigh(symm_pr)
            #assert scipy.allclose(scipy.linalg.inv(evecs), evecs.transpose())
            #assert scipy.allclose(symm_pr, scipy.dot(evecs, scipy.dot(scipy.diag(evals), evecs.transpose())))
            self.D[r] = evals
            self.Ainv[r] = evecs.transpose() * pr_half
            self.A[r] = (pr_neghalf * evecs.transpose()).transpose()

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
            _fill_diagonals(self.dPrxy['kappa'], self._diag_indices)
        if 'omega' in self.freeparams:
            scipy.copyto(self.dPrxy['omega'], self.Prxy / self.omega,
                    where=CODON_NONSYN)
            _fill_diagonals(self.dPrxy['omega'], self._diag_indices)
        if 'beta' in self.freeparams:
            self.dPrxy['beta'].fill(0)
            with scipy.errstate(divide='raise', under='raise', over='raise',
                    invalid='ignore'):
                scipy.copyto(self.dPrxy['beta'], self.Prxy *
                        (1 / self.beta + (self.piAx_piAy_beta *
                        scipy.log(self.piAx_piAy) / (1 - self.piAx_piAy_beta))),
                        where=CODON_NONSYN)
            scipy.copyto(self.dPrxy['beta'], self.Prxy/self.beta *
                    (1 - self.piAx_piAy_beta), where=(scipy.logical_and(
                    CODON_NONSYN, scipy.fabs(1 - self.piAx_piAy_beta)
                    < ALMOST_ZERO)))
            _fill_diagonals(self.dPrxy['beta'], self._diag_indices)
        if 'eta' in self.freeparams:
            for i in range(N_NT - 1):
                for w in range(i, N_NT):
                    scipy.copyto(self.dPrxy['eta'][i], self.Prxy / (self.eta[i]
                            - int(i == w)), where=CODON_NT_MUT[w])
                _fill_diagonals(self.dPrxy['eta'][i], self._diag_indices)

    def _update_B(self):
        """Update `B`."""
        for param in self.freeparams:
            if param == 'mu':
                continue
            paramval = getattr(self, param)
            if isinstance(paramval, float):
                self.B[param] = broadcastMatrixMultiply(self.Ainv,
                        broadcastMatrixMultiply(self.dPrxy[param], self.A))
            else:
                assert isinstance(paramval, numpy.ndarray) and paramval.ndim == 1
                for j in range(paramval.shape[0]):
                    self.B[param][j] = broadcastMatrixMultiply(self.Ainv,
                            broadcastMatrixMultiply(self.dPrxy[param][j],
                            self.A))

    def _update_dprx(self):
        """Update `dprx`."""
        if 'beta' in self.freeparams:
            for r in range(self.nsites):
                self.dprx['beta'][r] = self.prx[r] * (self.ln_pi_codon[r]
                        - scipy.dot(self.ln_pi_codon[r], self.prx[r]))
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
                        self.dprx['eta'][i][r] = self.prx[r] * (boolterm -
                                scipy.dot(boolterm, self.prx[r]) / self.prx[r].sum())

    def _fill_diagonals(self, m):
        """Fills diagonals of `nsites` matrices in `m` so rows sum to 0."""
        assert m.shape == (self.nsites, N_CODON, N_CODON)
        for r in range(self.nsites):
            scipy.fill_diagonal(m[r], 0)
            m[r][self._diag_indices] -= scipy.sum(m[r], axis=1)


class ExpCM_empirical_phi(ExpCM):
    """An `ExpCM` with `phi` calculated empirically from nucleotide frequencies.

    See `__init__` method for how to initialize an `ExpCM_empirical_phi`.

    The difference between and `ExpCM` and an `ExpCM_empirical_phi` is that
    the nucleotide frequency parameters `phi` are now calculated empirically
    from their observed frequences (denoted `g`) in the codon alignment such
    that the stationary state gives the empirically observed nucleotide
    frequencies.

    So compared to an `ExpCM`, `eta` (the transformation of `phi`) is no
    longer a free parameter, but is computed as a function of other
    model attributes.

    Has all the attributes of an `ExpCM` plus the following:
        `g` (`numpy.ndarray` of float, length `N_NT`)
            Empirical nucleotide frequencies in alignment, with
            `g[m]` being frequency of nucleotide `m`. Must sum
            to one.
        `dphi_dbeta` (`numpy.ndarray` of float, length `N_NT`)
            Derivative of `phi` with respect to `beta`.
    """

    # class variables
    _ALLOWEDPARAMS = copy.deepcopy(ExpCM._ALLOWEDPARAMS)
    _ALLOWEDPARAMS.remove('eta')
    _PARAMLIMITS = copy.deepcopy(ExpCM._PARAMLIMITS)
    _PARAMLIMITS['g'] = (0.05, 0.85)
    _PARAMTYPES = copy.deepcopy(ExpCM._PARAMTYPES)
    _PARAMTYPES['g'] = (scipy.ndarray, (N_NT,))

    def __init__(self, prefs, g, kappa=2.0, omega=0.5, beta=1.0, mu=1.0,
            freeparams=['kappa', 'omega', 'beta', 'mu']):
        """Initialize an `ExpCM_empirical_phi` object.

        Args:
            `prefs`, `kappa`, `omega`, `beta`, `mu`, `freeparams`
                Same meaning as for an `ExpCM`
            `g`
                Has the meaning described in the main class doc string.
        """

        _checkParam('g', g, self.PARAMLIMITS, self._PARAMTYPES)
        assert abs(1 - g.sum()) <= ALMOST_ZERO, "g doesn't sum to 1"
        self.g = g.copy()
        self.g /= self.g.sum()


        super(ExpCM_empirical_phi, self).__init__(prefs, kappa=kappa,
                omega=omega, beta=beta, mu=mu, freeparams=freeparams)

    def _update_phi(self):
        """Compute `phi`, `dphi_dbeta`, and `eta` from `g` and `frxy`."""
        self.phi = self._compute_empirical_phi(self.beta)
        _checkParam('phi', self.phi, self.PARAMLIMITS, self._PARAMTYPES)
        self._eta_from_phi()
        dbeta = 1.0e-3
        self.dphi_dbeta = scipy.misc.derivative(self._compute_empirical_phi,
                self.beta, dx=dbeta, n=1, order=5)
        dphi_dbeta_halfdx = scipy.misc.derivative(self._compute_empirical_phi,
                self.beta, dx=dbeta / 2, n=1, order=5)
        assert scipy.allclose(self.dphi_dbeta, dphi_dbeta_halfdx, atol=1e-5,
                rtol=1e-4), ("The numerical derivative dphi_dbeta differs "
                "considerably in value for step dbeta = {0} and a step "
                "half that size, giving values of {1} and {2}.").format(
                dbeta, self.dphi_dbeta, dphi_dbeta_halfdx)

    def _compute_empirical_phi(self, beta):
        """Returns empirical `phi` at the given value of `beta`.

        Does **not** set `phi` attribute, simply returns what
        should be value of `phi` given the current `g` and
        `pi_codon` attributes, plus the passed value of `beta`.
        Note that it uses the passed value of `beta`, **not**
        the current `beta` attribute.

        Initial guess is current value of `phi` attribute."""

        def F(phishort):
            """Difference between `g` and expected `g` given `phishort`."""
            phifull = scipy.append(phishort, 1 - phishort.sum())
            phiprod = scipy.ones(N_CODON, dtype='float')
            for w in range(N_NT):
                phiprod *= phifull[w]**CODON_NT_COUNT[w]
            frx_phiprod = frx * phiprod
            frx_phiprod_codonsum = frx_phiprod.sum(axis=1)
            gexpect = []
            for w in range(N_NT - 1):
                gexpect.append(
                        ((CODON_NT_COUNT[w] * frx_phiprod).sum(axis=1) /
                        frx_phiprod_codonsum).sum() / (3 * self.nsites))
            gexpect = scipy.array(gexpect, dtype='float')
            return self.g[ : -1] - gexpect

        frx = self.pi_codon**beta
        with scipy.errstate(invalid='ignore'):
            result = scipy.optimize.root(F, self.phi[ : -1].copy(),
                    tol=1e-8)
            assert result.success, "Failed: {0}".format(result)
            phishort = result.x
        return scipy.append(phishort, 1 - phishort.sum())

    def _update_dPrxy(self):
        """Update `dPrxy`, accounting for dependence of `phi` on `beta`."""
        super(ExpCM_empirical_phi, self)._update_dPrxy()
        if 'beta' in self.freeparams:
            self.dQxy_dbeta = scipy.zeros((N_CODON, N_CODON), dtype='float')
            for w in range(N_NT):
                scipy.copyto(self.dQxy_dbeta, self.dphi_dbeta[w],
                        where=CODON_NT_MUT[w])
            self.dQxy_dbeta[CODON_TRANSITION] *= self.kappa
            self.dPrxy['beta'] += self.Frxy * self.dQxy_dbeta
            _fill_diagonals(self.dPrxy['beta'], self._diag_indices)

    def _update_dprx(self):
        """Update `dprx`, accounting for dependence of `phi` on `beta`."""
        super(ExpCM_empirical_phi, self)._update_dprx()
        if 'beta' in self.freeparams:
            dphi_over_phi = scipy.zeros(N_CODON, dtype='float')
            for j in range(3):
                dphi_over_phi += (self.dphi_dbeta / self.phi)[CODON_NT_INDEX[j]]
            for r in range(self.nsites):
                self.dprx['beta'][r] += self.prx[r] * (dphi_over_phi
                        - scipy.dot(dphi_over_phi, self.prx[r]))


class ExpCM_empirical_phi_divpressure(ExpCM_empirical_phi):
    """`ExpCM` with empirical `phi` and *a priori* diversifying pressure at sites.

    See `__init__` method for how to initialize.

    Difference between `ExpCM_empirical_phi` and `ExpCM_empirical_phi_divpressure`
    is that `omega` is replaced by `omega * (1 + omega2 * deltar)` where
    `deltar` is the expectation of diversifying pressure at site `r`.

    The rate of non-synonymous to synonymous change is proportional to the
    expectation of diversifying pressure.

    Has all the attributes of an `ExpCM_empirical_phi` plus the following:
        `omega2` (float)
            Part of the expression which replaces `omega`.
        `divPressureValues` (`scipy.array` of length `nsites`)
            Array of `deltar` values in site order.
            The maximum absolute value of the array should be between 0 and 1.
    """

    # class variables
    _ALLOWEDPARAMS = copy.deepcopy(ExpCM_empirical_phi._ALLOWEDPARAMS)
    _ALLOWEDPARAMS.append('omega2')
    _REPORTPARAMS = copy.deepcopy(ExpCM_empirical_phi._REPORTPARAMS)
    _REPORTPARAMS.append('omega2')
    _PARAMLIMITS = copy.deepcopy(ExpCM_empirical_phi._PARAMLIMITS)
    _PARAMLIMITS['omega2'] = (-1, 999)
    _PARAMTYPES = copy.deepcopy(ExpCM_empirical_phi._PARAMTYPES)
    _PARAMTYPES['omega2'] = float

    def __init__(self, prefs, g, divPressureValues, kappa=2.0, omega=0.5,
            beta=1.0, mu=1.0,omega2=0.0,
            freeparams=['kappa', 'omega', 'beta', 'mu', 'omega2']):
        """Initialize an `ExpCM_empirical_phi_divpressure` object.

        Args:
            `prefs`, `kappa`, `omega`, `beta`, `mu`, `g`, `freeparams`
                Same meaning as for an `ExpCM_empirical_phi`
            `divPressureValues`, `omega2`
                Meaning described in the main class doc string.
        """
        _checkParam('omega2',omega2, self.PARAMLIMITS, self._PARAMTYPES)
        self.omega2 = omega2
        self.deltar = scipy.array(divPressureValues.copy())
        assert (max(scipy.absolute(self.deltar))) <= 1, (
                "A scaled deltar value is > 1 or < -1.")
        super(ExpCM_empirical_phi_divpressure, self).__init__(prefs, g,
                kappa=kappa, omega=omega, beta=beta, mu=mu,
                freeparams=freeparams)


    def _update_dPrxy(self):
        """Update `dPrxy`, accounting for dependence of `Prxy` on `omega2`."""
        super(ExpCM_empirical_phi_divpressure, self)._update_dPrxy()
        if 'omega2' in self.freeparams:
            with scipy.errstate(divide='raise', under='raise', over='raise',
                            invalid='ignore'):
                scipy.copyto(self.dPrxy['omega2'], -scipy.log(self.piAx_piAy_beta)
                        * self.Qxy * self.omega /
                        (1 - self.piAx_piAy_beta), where=CODON_NONSYN)
            scipy.copyto(self.dPrxy['omega2'], self.Qxy * self.omega,
                       where=(scipy.logical_and(CODON_NONSYN, scipy.fabs(1 -
                       self.piAx_piAy_beta) < ALMOST_ZERO)))
            for r in range(self.nsites):
                self.dPrxy['omega2'][r] *= self.deltar[r]
            _fill_diagonals(self.dPrxy['omega2'], self._diag_indices)


    def _update_Frxy(self):
        """Update `Frxy` from `piAx_piAy_beta`, `omega`, `omega2`, and `beta`."""
        with scipy.errstate(divide='raise', under='raise', over='raise',
                invalid='ignore'):
            scipy.copyto(self.Frxy, -1 * scipy.log(self.piAx_piAy_beta)
                    / (1 - self.piAx_piAy_beta), where=CODON_NONSYN)
        scipy.copyto(self.Frxy, 1, where=(scipy.logical_and(CODON_NONSYN,
                scipy.fabs(1 - self.piAx_piAy_beta) < ALMOST_ZERO)))
        for r in range(self.nsites):
            scipy.copyto(self.Frxy[r], self.Frxy[r] * self.omega *
                    (1 + self.omega2 * self.deltar[r]), where=CODON_NONSYN)

class YNGKP_M0(Model):
    """YNGKP_M0 model from Yang et al, 2000.

    See `__init__` method for how to initialize an `YNGKP_M0`.

    Attributes should **only** be updated via the `updateParams`
    method, do **not** set attributes directly.

    Attributes (see also those inherited from `Model`):
        `nsites` (int > 0)
            Number of amino acid sites in the gene
        `kappa` (float > 0)
            Transition-transversion ratio.
        `omega` (float > 0)
            Nonsynonymous to synonymous substitution ratio.
        `Pxy` (`numpy.ndarray` of floats, shape `(1, CODON, N_CODON)`
            `Pxy[0][x][y]` is substitution rate from codon `x` to `y`.
            Diagonal elements make rows sum to zero.
        `dPxy` (dict)
            Keyed by each string in `freeparams`, each value is `numpy.ndarray`
            of floats giving derivative of `Pxy` with respect to that parameter.
            The shape of each array `(1, N_CODON, N_CODON)`.
        `e_pw` (scipy.ndarray, shape `(3, N_NT)`)
            The empirical nucleotide frequencies for each position in a codon
            measured from the alignment. `e_pw[p][w]` give the frequency of
            nucleotide `w` at codon position `p`.
        `Phi_x` (scipy.ndarray, shape `(N_CODON,)`)
            The codon frequencies calculated from the `phi` frequencies.
            `Phi_x[x]` = phi[0][x0] * phi[1][x1] * phi[2][x2]
        `phi` (scipy.ndarray, shape `(3, N_NT)`)
            The model nucleotide frequencies for each position in a codon
            computed using the `CF3X4` estimator. `phi[p][w]` is the model
            frequency of nucleotide `w` at codon position `p`.

    Unlike the `ExpCM`, `YNGKP_M0` does not have site-specific calculations.
    In order to maintain consistency, most attributes, such as `Pxy` and `dPxy`
    do have a single "site" dimension which is carried through the calculations.
    When `M` and `dM` are returned, the single site matrix is repeated so the
    final dimensions of `M` and `dM` are match those in the docs for the `Models`
    abstract base class.

    The model nucleotide frequences, `phi` are computed using the
    empirical nucleotide frequencies `e_pw` and the `Corrected F3X4`
    estimator (Pond et al, 2010). Unlike the traditional `MG3X4` or `GY3X4`,
    this estimator takes into account stop codon nucleotide frequencies.
    The codon frequencies, `Phi_x` are computed from the `phi` values.
    """

    # class variables
    _REPORTPARAMS = ['kappa', 'omega', 'phi']
    _ALLOWEDPARAMS = ['kappa', 'omega', 'mu']
    _PARAMLIMITS = {'kappa':(0.01, 100.0),
                   'omega':(0.01, 100.0),
                   'mu':(1.0e-3, 1.0e3),
                   'e_pw':(0.02, 0.94),
                  }
    _PARAMTYPES = {'kappa':float,
                   'omega':float,
                   'mu':float,
                   'e_pw':(scipy.ndarray, (3, N_NT)),
                  }

    @property
    def ALLOWEDPARAMS(self):
        """See docs for `Model` abstract base class."""
        return self._ALLOWEDPARAMS

    @property
    def PARAMLIMITS(self):
        """See docs for `Model` abstract base class."""
        return self._PARAMLIMITS

    def __init__(self, e_pw, nsites, kappa=2.0, omega=0.5, mu=1.0,
            freeparams=['kappa', 'omega', 'mu']):
        """Initialize an `YNGKP_M0` object.

        Args:
            `kappa`, `omega`, `mu`,
                Model params described in main class doc string.
            `freeparams` (list of strings)
                Specifies free parameters.
            `e_pw`, `nsites`
                Meaning described in the main class doc string.
        """
        _checkParam('e_pw', e_pw, self.PARAMLIMITS, self._PARAMTYPES)
        self.e_pw = e_pw.copy()
        self.phi = self._calculate_correctedF3X4()
        assert scipy.allclose(self.phi.sum(axis = 1),\
                scipy.ones(3, dtype='float'),atol=1e-4, rtol=5e-3),\
                "The `phi` values do not sum to 1 for all `p`"
        # Not certain that this test should always be true,
        # so commenting it out.
        #assert (((self.e_pw - self.phi) * STOP_POSITIONS) >= 0.0).all(), (
        #        "Illogical `phi` values:\nphi = {0}\ne_pw = {1}\n"
        #        "(e_pw - phi) * STOP_POSITIONS = {2}".format(self.phi,
        #        self.e_pw, (self.e_pw - self.phi) * STOP_POSITIONS))

        self.Phi_x = scipy.ones(N_CODON, dtype='float')
        self._calculate_Phi_x()
        self._nsites = nsites
        assert self._nsites > 0, "There must be more than 1 site in the gene"

        #check allowed params
        assert all(map(lambda x: x in self.ALLOWEDPARAMS, freeparams)),\
                "Invalid entry in freeparams\nGot: {0}\nAllowed: {1}".format(
                ', '.join(freeparams), ', '.join(self.ALLOWEDPARAMS))
        self._freeparams = list(freeparams) # underscore as `freeparams` is property

        # set attributes to calling params
        self._mu = mu # underscore as `mu` is property
        self.kappa = kappa
        self.omega = omega
        for (name, value) in [('kappa', self.kappa), ('omega', self.omega),
                    ('mu', self.mu)]:
            _checkParam(name, value, self.PARAMLIMITS, self._PARAMTYPES)

        # define other params, initialized appropriately
        #single site dimension to be carried through the calcs added here
        self.Pxy = scipy.zeros((1, N_CODON, N_CODON), dtype='float')
        self.D = scipy.zeros((1, N_CODON), dtype='float')
        self.A = scipy.zeros((1, N_CODON, N_CODON), dtype='float')
        self.Ainv = scipy.zeros((1, N_CODON, N_CODON), dtype='float')
        self.dPxy = {}
        self.B = {}
        for param in self.freeparams:
            if param in self._ALLOWEDPARAMS:
                self.dPxy[param] = scipy.zeros((1, N_CODON, N_CODON),
                        dtype='float')
                self.B[param] = scipy.zeros((1, N_CODON, N_CODON),
                        dtype='float')
            else:
                raise ValueError("Unrecognized param {0}".format(param))

        # indexes diagonals in square matrices
        self._diag_indices = scipy.diag_indices(N_CODON)
        self.updateParams({}, update_all=True)

    @property
    def stationarystate(self):
        """See docs for `Model` abstract base class."""
        #repeat the single state by the number sites
        return scipy.tile(self.Phi_x, (self.nsites, 1)) 

    def dstationarystate(self, param):
        """See docs for `Model` abstract base class."""
        return scipy.zeros((self.nsites, N_CODON), dtype='float')

    @property
    def paramsReport(self):
        """See docs for `Model` abstract base class."""
        report = {}
        for param in self._REPORTPARAMS:
            pvalue = getattr(self, param)
            if isinstance(pvalue, float):
                report[param] = pvalue
            elif isinstance(pvalue, scipy.ndarray) and pvalue.shape == (3, N_NT):
                for p in range(3):
                    for w in range(N_NT - 1):
                        report['{0}{1}{2}'.format(param, p, INDEX_TO_NT[w])] =\
                                pvalue[p][w]
            else:
                raise RuntimeError("Unexpected param: {0}".format(param))
        return report

    @property
    def branchScale(self):
        """See docs for `Model` abstract base class."""
        bs = -(self.Phi_x * scipy.diagonal(self.Pxy[0])).sum() * self.mu
        assert bs > 0
        return bs

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

    @mu.setter
    def mu(self, value):
        """Set new `mu` value."""
        self._mu = value

    def _calculate_correctedF3X4(self):
        '''Calculate `phi` based on the empirical `e_pw` values'''
        def F(phi):
            phi_reshape = phi.reshape((3, N_NT))
            functionList = []
            stop_frequency = []

            for x in range(N_STOP):
                codonFrequency = STOP_CODON_TO_NT_INDICES[x] * phi_reshape
                codonFrequency = scipy.prod(codonFrequency.sum(axis=1))
                stop_frequency.append(codonFrequency)
            C = scipy.sum(stop_frequency)

            for p in range(3):
                for w in range(N_NT):
                    s = 0
                    for x in range(N_STOP):
                        if STOP_CODON_TO_NT_INDICES[x][p][w] == 1:
                            s += stop_frequency[x]
                    functionList.append((phi_reshape[p][w] - s)/(1 - C)
                            - self.e_pw[p][w])
            return functionList

        phi = self.e_pw.copy().flatten()
        with scipy.errstate(invalid='ignore'):
            result = scipy.optimize.root(F, phi,
                    tol=1e-8)
            assert result.success, "Failed: {0}".format(result)
            return result.x.reshape((3, N_NT))

    def _calculate_Phi_x(self):
        """Calculate `Phi_x` (stationary state) from `phi`."""
        self.Phi_x = scipy.ones(N_CODON, dtype='float')
        for codon in range(N_CODON):
            for pos in range(3):
                self.Phi_x[codon] *= self.phi[pos][CODON_NT_INDEX[pos][codon]]

    def updateParams(self, newvalues, update_all=False):
        """See docs for `Model` abstract base class."""
        assert all(map(lambda x: x in self.freeparams, newvalues.keys())),\
                "Invalid entry in newvalues: {0}\nfreeparams: {1}".format(
                ', '.join(newvalues.keys()), ', '.join(self.freeparams))
        changed = set([]) # contains string names of changed params
        for (name, value) in newvalues.items():
            _checkParam(name, value, self.PARAMLIMITS, self._PARAMTYPES)
            if isinstance(value, scipy.ndarray):
                if (value != getattr(self, name)).any():
                    changed.add(name)
                    setattr(self, name, value)
            else:
                if value != getattr(self, name):
                    changed.add(name)
                    setattr(self, name, value)

        if update_all or changed:
            self._cached = {}

        # The order of the updating below is important.
        # If you change it, you may break either this class
        # **or** classes that inherit from it.
        # Note also that not all attributes need to be updated
        # for all possible parameter changes, but just doing it
        # this way is much simpler and adds negligible cost.
        if update_all or (changed and changed != set(['mu'])):
            self._update_Pxy()
            self._update_Pxy_diag()
            self._update_dPxy()
            self._update_B()

    def M(self, t, tips=None, gaps=None):
        """See docs for method in `Model` abstract base class."""
        assert isinstance(t, float) and t > 0, "Invalid t: {0}".format(t)
        with scipy.errstate(under='ignore'): # don't worry if some values 0
            if ('expD', t) not in self._cached:
                self._cached[('expD', t)] = scipy.exp(self.D * self.mu * t)
            expD = self._cached[('expD', t)]
            # swap axes to broadcast multiply D as diagonal matrix
            temp = scipy.ascontiguousarray((self.A.swapaxes(0, 1) * expD).swapaxes(1, 0), dtype = float)
            M = broadcastMatrixMultiply(temp, self.Ainv)
            if tips is None:
                return scipy.tile(M, (self.nsites,1,1))
            else:
                newM = scipy.zeros((len(tips), N_CODON))
                for i in range(len(tips)):
                    newM[i] =(M[0][:,tips[i]])
                if gaps is not None:
                    newM[gaps] = scipy.ones(N_CODON, dtype='float')
                return newM

    def dM(self, t, param, Mt, tips=None, gaps=None):
        """See docs for method in `Model` abstract base class."""
        assert isinstance(t, float) and t > 0, "Invalid t: {0}".format(t)
        assert param in self.freeparams, "Invalid param: {0}".format(param)
        if param == 'mu':
            if tips is None:
                dM_param = scipy.tile(broadcastMatrixMultiply(self.Pxy, scipy.tile(Mt[0], (1,1,1)), alpha=t), (self.nsites, 1, 1))
            else:
                #Pxy is tiled over the number of sites
                dM_param = broadcastMatrixVectorMultiply(scipy.tile(self.Pxy[0], (self.nsites,1,1)), Mt, alpha=t)
                if gaps is not None:
                    dM_param[gaps] = scipy.zeros(N_CODON, dtype='float')
            return dM_param

        paramval = getattr(self, param)
        assert isinstance(paramval, float), "All params should be floats"

        if ('expD', t) not in self._cached:
            self._cached[('expD', t)] = scipy.exp(self.D * self.mu * t)
        expD = self._cached[('expD', t)]

        if ('V', t) not in self._cached:
            if 'Dxx_Dyy' not in self._cached:
                Dyy = scipy.tile(self.D, (1, N_CODON)).reshape(
                        1, N_CODON, N_CODON)
                Dxx = scipy.array([Dyy[r].transpose() for r in
                        range(1)])
                self._cached['Dxx_Dyy'] = Dxx - Dyy
            Dxx_Dyy = self._cached['Dxx_Dyy']
            if 'Dxx_Dyy_lt_ALMOST_ZERO' not in self._cached:
                self._cached['Dxx_Dyy_lt_ALMOST_ZERO'] = scipy.fabs(
                        Dxx_Dyy) < ALMOST_ZERO
            Dxx_Dyy_lt_ALMOST_ZERO = self._cached['Dxx_Dyy_lt_ALMOST_ZERO']
            with scipy.errstate(divide='raise', under='ignore',
                    over='raise', invalid='ignore'):
                expDyy = scipy.tile(expD,(1, N_CODON)).reshape(
                        1, N_CODON, N_CODON)
                expDxx = scipy.array([expDyy[r].transpose() for r in
                        range(1)])
                V = (expDxx - expDyy) / Dxx_Dyy
            with scipy.errstate(under='ignore'): # OK if some values 0
                scipy.copyto(V, self.mu * t * expDxx, where=
                        Dxx_Dyy_lt_ALMOST_ZERO)
            self._cached[('V', t)] = V
        V = self._cached[('V', t)]

        with scipy.errstate(under='ignore'): # don't worry if some values 0
            dM_param = broadcastMatrixMultiply(self.A,
                        broadcastMatrixMultiply(self.B[param]
                        * V, self.Ainv))
            if tips is None:
                return scipy.tile(dM_param, (self.nsites, 1, 1))
            else:
                newdM_param = scipy.zeros((len(tips), N_CODON))
                for i in range(len(tips)):
                    newdM_param[i] =(dM_param[0][:,tips[i]])
                if gaps is not None:
                    newdM_param[gaps] = scipy.zeros(N_CODON, dtype='float')
                return newdM_param

    def _update_Pxy(self):
        """Update `Pxy` using current `omega`, `kappa`, and `Phi_x`."""
        scipy.copyto(self.Pxy, self.Phi_x.transpose(), where=CODON_SINGLEMUT)
        self.Pxy[0][CODON_NONSYN] *= self.omega
        self.Pxy[0][CODON_TRANSITION] *= self.kappa
        _fill_diagonals(self.Pxy, self._diag_indices)

    def _update_dPxy(self):
        """Update `dPxy`."""
        if 'kappa' in self.freeparams:
            scipy.copyto(self.dPxy['kappa'], self.Pxy / self.kappa,
                    where=CODON_TRANSITION)
            _fill_diagonals(self.dPxy['kappa'], self._diag_indices)
        if 'omega' in self.freeparams:
            scipy.copyto(self.dPxy['omega'], self.Pxy / self.omega,
                    where=CODON_NONSYN)
            _fill_diagonals(self.dPxy['omega'], self._diag_indices)

    def _update_Pxy_diag(self):
        """Update `D`, `A`, `Ainv` from `Pxy`, `Phi_x`."""
        for r in range(1):
            Phi_x_half = self.Phi_x**0.5
            Phi_x_neghalf = self.Phi_x**-0.5
            #symm_p = scipy.dot(scipy.diag(Phi_x_half), scipy.dot(self.Pxy[r], scipy.diag(Phi_x_neghalf)))
            symm_p = (Phi_x_half * (self.Pxy[r] * Phi_x_neghalf).transpose()).transpose()
            #assert scipy.allclose(symm_p, symm_p.transpose())
            (evals, evecs) = scipy.linalg.eigh(symm_p)
            #assert scipy.allclose(scipy.linalg.inv(evecs), evecs.transpose())
            #assert scipy.allclose(symm_pr, scipy.dot(evecs, scipy.dot(scipy.diag(evals), evecs.transpose())))
            self.D[r] = evals
            self.Ainv[r] = evecs.transpose() * Phi_x_half
            self.A[r] = (Phi_x_neghalf * evecs.transpose()).transpose()

    def _update_B(self):
        """Update `B`."""
        for param in self.freeparams:
            if param == 'mu':
                continue
            paramval = getattr(self, param)
            assert isinstance(paramval, float), "Paramvalues must be floats"
            self.B[param] = broadcastMatrixMultiply(self.Ainv,
                    broadcastMatrixMultiply(self.dPxy[param], self.A))


def _checkParam(param, value, paramlimits, paramtypes):
    """Checks if `value` is allowable value for `param`.

    Raises except if `value` is not acceptable, otherwise
    returns `None` if value is acceptable.

    `paramlimits` and `paramtypes` are the `PARAMLIMITS`
    and `_PARAMTYPES` attributes of a `Model`.
    """
    assert param in paramlimits, "Invalid param: {0}".format(param)
    (lowlim, highlim) = paramlimits[param]
    paramtype = paramtypes[param]
    if isinstance(paramtype, tuple):
        (paramtype, paramshape) = paramtype
        if not (isinstance(value, paramtype)):
            raise ValueError("{0} must be {1}, not {2}".format(
                    param, paramtype, type(param)))
        if value.shape != paramshape:
            raise ValueError("{0} must have shape {1}, not {2}".format(
                    param, paramshape, value.shape))
        if value.dtype != 'float':
            raise ValueError("{0} must have dtype float, not {1}".format(
                    param, value.dtype))
        if not ((lowlim <= value).all() and (value <= highlim).all()):
            raise ValueError("{0} must be >= {1} and <= {2}, not {3}".format(
                    param, lowlim, highlim, value))
    else:
        if not isinstance(value, paramtype):
            raise ValueError("{0} must be a {1}, not a {2}".format(
                    param, paramtype, type(value)))
        if not (lowlim <= value <= highlim):
            raise ValueError("{0} must be >= {1} and <= {2}, not {3}".format(
                    param, lowlim, highlim, value))


def _fill_diagonals(m, diag_indices):
    """Fills diagonals of `nsites` matrices in `m` so rows sum to 0."""
    assert m.ndim ==3, "M must have 3 dimensions"
    assert m.shape[1] == m.shape[2], "M must contain square matrices"
    for r in range(m.shape[0]):
        scipy.fill_diagonal(m[r], 0)
        m[r][diag_indices] -= scipy.sum(m[r], axis=1)


def DiscreteGamma(alpha, beta, ncats):
    """Returns category means for discretized gamma distribution.

    The distribution is evenly divided into categories, and the
    mean of each category is returned.

    Args:
        `alpha` (`float` > 0)
            Shape parameter of gamma distribution.
        `beta` (`float` > 0)
            Inverse scale parameter of gamma distribution.
        `ncats` (`int` > 0)
            Number of categories in discretized gamma distribution.

    Returns:
        `catmeans` (`scipy.ndarray` of `float`, shape `(ncats,)`)
            `catmeans[k]` is the mean of category `k` where 
            `0 <= k < ncats`.

    Check that we get values in Figure 1 of Yang, J Mol Evol, 39:306-314
    >>> catmeans = DiscreteGamma(0.5, 0.5, 4)
    >>> scipy.allclose(catmeans, scipy.array([0.0334, 0.2519, 0.8203, 2.8944]), atol=1e-4)
    True

    Make sure we get expected mean of alpha / beta
    >>> alpha = 0.6
    >>> beta = 0.8
    >>> ncats = 6
    >>> catmeans = DiscreteGamma(alpha, beta, ncats)
    >>> scipy.allclose(catmeans.sum() / ncats, alpha / beta)
    True
    """
    assert alpha > 0
    assert beta > 0
    assert ncats > 0
    scale = 1.0 / beta
    catmeans = scipy.ndarray(ncats, dtype='float')
    for k in range(ncats):
        if k == 0:
            lower = 0.0
            gammainc_lower = 0.0
        else:
            lower = upper
            gammainc_lower = gammainc_upper
        if k == ncats - 1:
            upper = float('inf')
            gammainc_upper = 1.0
        else:
            upper = scipy.stats.gamma.ppf((k + 1) / float(ncats), alpha, 
                    scale=scale)
            gammainc_upper = scipy.special.gammainc(alpha + 1, upper * beta)
        catmeans[k] = alpha * ncats * (gammainc_upper - gammainc_lower) / beta
    assert scipy.allclose(catmeans.sum() / ncats, alpha / beta), (
            "catmeans is {0}, mean of catmeans is {1}, expected mean "
            "alpha / beta = {2} / {3} = {4}").format(catmeans, 
            catmeans.sum() / ncats, alpha, beta, alpha / beta)
    return catmeans


if __name__ == '__main__':
    import doctest
    doctest.testmod()
