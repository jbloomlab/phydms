"""Module for computing difference in preferences by optimizing preferences with prior."""


import sys
import math
import numpy
import scipy.optimize
import scipy.stats



class PrefsToVec(object):
    """Converts dictionary of preferences to vector of values between 0 and 1.

    The preferences are most easily stored as a dictionary keyed by character
    with values that sum to one. However, for optimization it is helpful to
    re-represent these as a vector of numbers constrained to have values
    between 0 and 1 (inclusive). 

    Specifically, let the preferences for the different sites by
    :math:`\pi_i` where :math:`1 = \sum_i \pi_i` and 
    :math:`0 \le i \lt n`. Let :math:`x_i` be the vector element
    for :math:`i \lt n - 1`. We define 
    :math:`\pi_i = x_i \prod_{j \lt i} \left(1 - x_j\\right)`
    for :math:`i \lt n - 1`, and 
    :math:`\pi_{n-1} = 1 - \sum_{j = 0}^{n - 2} \pi_j`.
    This implies a definition for :math:`x_i` of
    :math:`x_i = \pi_i / \left(1 - \sum_{j \lt i} \pi_j\\right)`

    To set up an object, initialize with::

        prefstovec = PrefsToVec(initprefs)

    where *initprefs* is some initial set of preferences as a dictionary
    keyed by all characters. The value of *initprefs* is used to sort
    the indices in the vector (they always go from largest to smallest
    element, which might help avoid numerical underflow for last values).

    To subsequently convert dictionaries of preferences vectors, use::

        vec = prefstovec.Vector(prefs)

    and to convert back::

        prefs = prefstovec.Prefs(vec)

    Here is an example:

    >>> initprefs = {'A':0.45, 'C':0.25, 'G':0.21, 'T':0.09}
    >>> prefstovec = PrefsToVec(initprefs)
    >>> prefs = {'A':0.4, 'C':0.3, 'G':0.2, 'T':0.1}
    >>> vec = prefstovec.Vector(prefs)
    >>> numpy.allclose(vec, numpy.array([0.4, 0.5, 2. / 3.]))
    True
    >>> prefs2 = prefstovec.Prefs(vec)
    >>> all([abs(prefs[x] - prefs2[x]) < 1.e-5 for x in prefs.keys()])
    True
    """

    def __init__(self, initprefs, tol=1e-6):
        """Initializes class with *initprefs* as starting prefs dictionary.
        
        *tol* is the tolerance for preferences summing to one."""
        assert tol < 1e-4, "Unreasonably large tol: %s" % tol
        self.tol = tol
        # dicts to convert from prefs chars (e.g. amino acids) to vector indices
        # order from largest to smallest
        decorated_chars = [(ipi, ichar) for (ichar, ipi) in initprefs.items()]
        decorated_chars.sort()
        decorated_chars.reverse()
        self.chars = [tup[1] for tup in decorated_chars]
        self.nchars = len(self.chars)
        assert self.nchars > 1, "Need more than one character"

    def Vector(self, prefs):
        """Converts prefs dictionary to vector."""
        assert len(prefs) == self.nchars, "prefs not for right number of characters"
        assert all([pi >= 0 for pi in prefs.values()]), "prefs not all >= 0: %s" % str(prefs)
        assert abs(1.0 - sum(prefs.values())) < self.tol, "prefs do not sum to one: %s" % str(prefs)
        assert set(self.chars) == set(prefs.keys()), "prefs does not have correct character keys: %s" % str(prefs)
        vec = []
        runningsum = 0.0
        for char in self.chars[ : -1]:
            if runningsum >= 1.0 - self.tol:
                vec.append(0.5)
            else:
                vec.append(min(1.0, prefs[char] / (1.0 - runningsum)))
            runningsum += prefs[char]
        return numpy.array(vec)

    def Prefs(self, vec):
        """Converts vector to prefs dictionary.
        
        All values are adjusted to be at least equal to the *tol*
        used when initializing this object."""
        assert len(vec) == self.nchars - 1, "vec not for right number of characters"
        assert all(vec >= 0) and all(vec <= 1 + self.tol), "vec elements are not all >= 0 and <= 1"
        prefs = {}
        runningprod = 1.0
        for (ichar, char) in enumerate(self.chars[ : -1]):
            veci = min(1.0, vec[ichar])
            prefs[char] = veci * runningprod
            runningprod *= (1.0 - veci)
        prefs[self.chars[-1]] = max(0, 1.0 - sum(prefs.values()))
        maxprefchar = [(pi, char) for (char, pi) in prefs.items()]
        maxprefchar.sort()
        maxprefchar = maxprefchar[-1][1]
        for char in prefs.keys():
            if prefs[char] < self.tol:
                prefs[char] += self.tol
                prefs[maxprefchar] -= self.tol
        assert abs(1.0 - sum(prefs.values())) < self.tol, "prefs do not sum to one: %s\n%s" % (str(prefs), str(vec))
        assert all([pi >= self.tol for pi in prefs.values()]), "prefs has element less than tol = %g:\n%s" % (self.tol, str(prefs))
        return prefs


class InvQuadPrefsPrior(object):
    """Computes inverse quadratic prior over preferences.

    This object is instantiated like this::

        prefsprior = InvQuadPrefsPrior(peakprefs, cparams)

    where *peakprefs* is a dictionary of the preferences that specifies
    the peak of the prior, and *cparams* is the tuple
    *(c1, c2)* of numbers >= 0
    that specify how "concentrated" the prior is (larger value means
    more concentrated around *peakprefs*). 
    
    *tol* is the tolerance
    for having preferences not sum to one.

    The prior is defined as follows. Given some new vector with elements
    :math:`\hat{\pi}_{r,a}` and initial (*peakprefs*) values
    of :math:`\pi_{r,a}`, the log prior probability is
    :math:`\log\left[\Pr\left(\left\{\hat{\pi}_{r,a}\\right\} \mid \left\{\pi_{r,a}\\right\}, \\beta\\right)\\right] = - C_2 \sum_a \log\left(1 + C_1 \\times \left(\hat{\pi}_{r,a} - \left(\pi_{r,a}\\right)^{\\beta}\\right)^2\\right)`.

    To return the log prior probability of any given set of preferences
    contained in a dictionary *prefs*, simply use::

        prior = prefsprior.LogPrior(prefs)

    """
    def __init__(self, peakprefs, cparams, tol=1e-6):
        """Initialize object."""
        (self.c1, self.c2) = cparams
        assert self.c1 >= 0, "c1 must be >= 0"
        assert self.c2 >= 0, "c2 must be >= 0"
        assert tol < 1e4, "Unreasonably large tol: %s" % tol
        self.tol = tol
        assert abs(1.0 - sum(peakprefs.values())) < tol, "peakprefs do not sum to one:\n%s" % str(peakprefs)
        assert all([pi >= 0 for pi in peakprefs.values()]), "peakprefs not all >= 0: %s" % str(peakprefs)
        self.chars = peakprefs.keys()
        self.priorvec = numpy.array([peakprefs[char] for char in self.chars])

    def LogPrior(self, prefs):
        """Return prior probability of *prefs* (which is in dictionary format)."""
        assert set(prefs.keys()) == set(self.chars), "Invalid character keys in prefs:\n%s" % str(prefs)
        assert abs(1.0 - sum(prefs.values())) < self.tol, "prefs do not sum to one:\n%s" % str(prefs)
        assert all([pi >= 0 for pi in prefs.values()]), "prefs not all >= 0: %s" % str(prefs)
        prefsvec = numpy.array([prefs[char] for char in self.chars])
        logp = 0.0
        for (pi1, pi2) in zip(self.priorvec, prefsvec):
            logp -= self.c2 * math.log(1.0 + self.c1 * (pi1 - pi2)**2)
        return logp


class DirichletPrefsPrior(object):
    """Computes dirichlet prior over preferences.

    This object is instantiated like this::

        prefsprior = DirichletPrefsPrior(peakprefs, concentration, minvalue)

    where *peakprefs* is a dictionary of the preferences that specifies
    the peak of the prior, and *concentration* is a number > 1 that
    specifies how "concentrated" the prior is (larger value means
    more concentrated). *minvalue* gives the minimum allowed value
    of the peak preferences for any character; it should be a small number
    but > 0. Any preference in *peakprefs* that is less than *minvalue* 
    is adjusted up to be equal to *minvalue*. *tol* is the tolerance
    for having preferences not sum to one.

    The prior itself is a Dirichlet. It turns out that if a Dirichlet 
    is parameterized by a vector of values :math:`\\alpha_i`, then
    the mode for element :math:`i` is :math:`\\frac{\\alpha_i - 1}{\sum_j \\alpha_j - n}`
    where :math:`n` is the number of elements in the vector. Let :math:`C > 1`
    denote the concentration parameter *concentration*. We can parameterize
    a Dirichlet so that it has the mode :math:`\pi_i` by choosing a parameter
    vector with elements
    :math:`\\alpha_i = \left(C - 1\\right) \\times n \\times \pi_i + 1`.

    To return the log prior probability of any given set of preferences
    contained in a dictionary *prefs*, simply use::

        prior = prefsprior.LogPrior(prefs)

    Here is an example:

    >>> peakprefs = {'A':0.4, 'C':0.3, 'G':0.3, 'T':0.0}
    >>> concentration = 2
    >>> minvalue = 1.0e-3
    >>> prefsprior = DirichletPrefsPrior(peakprefs, concentration, minvalue)
    >>> prefs1 = {'A':0.39, 'C':0.3, 'G':0.3, 'T':0.01}
    >>> prefs2 = {'A':0.09, 'C':0.3, 'G':0.3, 'T':0.31}
    >>> prefsprior.LogPrior(prefs1) > prefsprior.LogPrior(prefs2)
    True
    """

    def __init__(self, peakprefs, concentration, minvalue, tol=1e-6):
        """Initialize object."""
        assert concentration > 1, "concentration should be > 1"
        assert 1e-2 > minvalue > 0, "Unreasonable value of %g for minvalue"
        assert tol < 1e4, "Unreasonably large tol: %s" % tol
        self.tol = tol
        assert abs(1.0 - sum(peakprefs.values())) < tol, "peakprefs do not sum to one:\n%s" % str(peakprefs)
        assert all([pi >= 0 for pi in peakprefs.values()]), "peakprefs not all >= 0: %s" % str(peakprefs)
        self.chars = peakprefs.keys()
        self.priorvec = numpy.array([(concentration - 1.0) * len(self.chars) * max(minvalue, peakprefs[char]) + 1.0 for char in self.chars])
        self.dirichlet = scipy.stats.dirichlet(self.priorvec)

    def LogPrior(self, prefs):
        """Return prior probability of *prefs* (which is in dictionary format)."""
        assert set(prefs.keys()) == set(self.chars), "Invalid character keys in prefs:\n%s" % str(prefs)
        assert abs(1.0 - sum(prefs.values())) < self.tol, "prefs do not sum to one:\n%s" % str(prefs)
        assert all([pi >= 0 for pi in prefs.values()]), "prefs not all >= 0: %s" % str(prefs)
        prefsvec = [prefs[char] for char in self.chars[ : -1]]
        return self.dirichlet.logpdf(prefsvec)



def OptimizePrefs(tl, prefs, site, prior, concentration, minvalue=1e-4, noprior=False, nologl=False, randinitprefs=False, optmethod='SLSQP'):
    """Optimize preferences along with a prior constraint.

    *tl* is a *phydmslib.pybpp.PyBppTreeLikelihood* object that has
    an experimentally defined codon model for site *site*. Generally,
    all parameters should alread be at their desired values
    except the preferences, which are sampled by MCMC.

    *prefs* gives the initial estimates for the preferences that is
    used to center the prior over the inferred preferences. 
    *prefs[x]* should gives the estimate for the preference for 
    character *x* at site *site* (*prefs* holds preferences just
    for this one site).

    *site* is an integer giving the site in *tl* for which we
    performing MCMC.

    *prior* is the type of prior used over the preferences. Must be
    one of the following string:

        - *invquad* : the prior described in the documentation
          for *InvQuadPrefsPrior*.

        - *dirichlet* : the prior described in the documentation 
          for *DirichletPrefsPrior*.

    *concentration* is how strongly we concentrate the prior.
    Larger values lead to priors more strongly peaked on the initial
    value of *site*. Should be a value appropriate for passing to 
    *InvQuadPrefsPrior* (if using *prior* of *invquad*) or to
    *DirichletPrefsPrior* (if using *prior* of *dirichlet*).
    
    *minvalue* is a lower bound on the minimum valued allowed for the
    peak prior estimate for any preference and for any preference. 
    is *dirichlet*.

    *noprior* and *nologl* are for debugging only; the indicate that we don't
    include the prior or the likelihood in the optimization.

    *randinitprefs* specify that we seed the optimization at random preferences
    chosen from the seed *randinitprefs*, which should be an integer. Otherwise
    if *False*, we seed optimization at *prefs*.

    *optmethod* is the *scipy.optimize.minimize* method to use. Should be one that
    accepts bounds (constrained minimization).

    The return value is *(optimizedprefs, optstring, converged)*.
    *optimizedprefs* is in the same format as *prefs* and gives
    the optimized preferences. *optstring* describes the convergence
    results, and *converged* is *True* if the optimization converged.
    """
    # error check arguments
    assert (not noprior) or (not nologl), "Cannot use both noprior and logl"
    assert all([pi >= 0 for pi in prefs.values()]), "The preference must be > 0 for all sites"
    prefs = dict([(aa, max(pi, 1.0e-7)) for (aa, pi) in prefs.items()]) # make all a bit greater than zero
    assert abs(sum(prefs.values()) - 1.0) < 1.0e-5, "The sum of the preferences must be one"
    assert isinstance(site, int) and 1 <= site <= tl.NSites(), "site of %d is not in the tree likelihood object" % site
    assert set(prefs.keys()) == set(tl.GetPreferences(site).keys()), "prefs does not have keys for the same characters as the tree likelihood object"
    assert 0 < minvalue < 1e-2, "Unreasonable value for minvalue: %g" % minvalue

    # set up initial preference vector, prior vector, and the bounds
    prefstovec = PrefsToVec(prefs)
    if randinitprefs != False:
        numpy.random.seed(randinitprefs)
        aas = prefs.keys()
        initprefs = dict(zip(aas, numpy.random.dirichlet([len(aas) * prefs[aa] for aa in aas])))
        initvec = prefstovec.Vector(initprefs)
    else:
        initvec = prefstovec.Vector(prefs)
    bounds = [(minvalue, 1.0 - minvalue) for i in range(len(initvec))]
    if prior == 'invquad':
        prefsprior = InvQuadPrefsPrior(prefs, concentration)
    elif prior == 'dirichlet':
        prefsprior = DirichletPrior(prefs, concentration, minvalue)
    else:
        raise ValueError("Invalid value of prior: %s" % prior)

    # function to minimize
    def NegLogPosterior(vec):
        """Returns **negative** log likelihood as we are minimizing."""
        if any(numpy.isnan(vec)):
            raise RuntimeError("Call outside support: %s" % str(vec))
        else:
            iprefs = prefstovec.Prefs(vec)
            if any([pi <= 0 for pi in iprefs.values()]) or any([pi >= 1 for pi in iprefs.values()]):
                raise RuntimeError("iprefs outside support:\n%s\n%s" % (str(iprefs), str(vec)))
            if noprior and nologl:
                raise ValueError("Can't use noprior and nologl")
            elif noprior:
                tl.SetPreferences(iprefs, site)
                return -(tl.LogLikelihood())
            elif nologl:
                return -(prefsprior.LogPrior(iprefs))
            else:
                tl.SetPreferences(iprefs, site)
                return -(tl.LogLikelihood() + prefsprior.LogPrior(iprefs))

    # do the minimization
    result = scipy.optimize.minimize(NegLogPosterior, initvec, method=optmethod, bounds=bounds)
    
    # process the results
    assert len(result.x) == len(initvec), "result has length that differs from initvec"
    optimizedprefs = prefstovec.Prefs(result.x)
    assert all([ipi >= 0 for ipi in optimizedprefs.values()]) and abs(1.0 - sum(optimizedprefs.values())) < 1.0e5, "Invalid optimized prefs: %s" % str(optimizedprefs)
    return (optimizedprefs, result.message, result.success)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
