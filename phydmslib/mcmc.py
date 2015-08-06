"""Module for MCMC; uses ``pymc`` (tested with version 2.3.4)."""


import math
import random
import pymc


def PrefsMCMC(tl, prefs, site, concentrationparam, seed, verbose=0, nchains=2, nsteps=1e4, burnfrac=0.2, nincreasetries=2, foldincrease=2, convergence=1.1):
    """MCMC to infer posterior mean preferences from tree likelihood.

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

    *concentrationparam* is the concentration parameter for the 
    Dirichlet prior centered at *prefs*. The element in the 
    Dirichlet vector for site *x* is the product of
    *concentrationparam*, *prefs[x]*, and *len(prefs)*. So for
    a uniform Dirichlet (all preferences equal), a value of
    *concentrationparam* equal to one corresponds to equal
    prior probability for all preference vectors.

    *seed* is the random number seed.

    *verbose* is the MCMC verbosity. Smaller integers mean less
    output (0 is no output), you will get more output with values
    up to at least 3.

    *nchains*, *burnfrac*, *nsteps*, *nincreasetries*, *foldincrease*, 
    *convergence* specify how the MCMC is performed. The MCMC
    is considered to converge if the chains (*nchains* must be >= 2)
    have a Gelman-Rubin R that is less than *convergence* **and**
    if the max root-mean-square difference in the inferred preferences
    between chains is less than *(convergence - 1.0) / 2*. First, we try
    *nchains* chains with *nsteps* steps, of which the first *burnfrac*
    is burn-in. If convergence fails, we try again *nincreasetries* more
    times, each time increasing the number of steps by *foldincrease*.
    Each chain is always seeded with preferences drawn
    randomly from a flat Dirichlet distribution.

    The return value is *(meanprefs, mcmcstring, converged)* where:

        - *meanprefs* is a dictionary with the same keys as *prefs*,
          and with values giving the posterior mean preference for
          that character.

        - *mcmcstring* is a string that summarizes some information
          about the MCMC and its convergence.

        - *converged* is *True* if the MCMC converged and *False* 
          otherwise.
    """ 
    # check variables
    assert nchains > 1
    assert convergence > 1.0
    assert foldincrease > 1
    assert nincreasetries >= 0
    assert concentrationparam > 0, "concentratioparam must be > 0"
    assert all([pi > 0 for pi in prefs.values()]), "The preference must be > 0 for all sites"
    assert abs(sum(prefs.values()) - 1.0) < 1.0e-5, "The sum of the preferences must be one"
    assert isinstance(site, int) and 1 <= site <= tl.NSites(), "site of %d is not in the tree likelihood object" % site
    assert set(prefs.keys()) == set(tl.GetPreferences(site).keys()), "prefs does not have keys for the same characters as the tree likelihood object"

    # seed random number generator
    pymc.numpy.random.seed(seed)

    # dicts that convert from prefs chars (e.g. amino acids) to vector indices and back
    chars = prefs.keys()
    chars.sort()
    nchars = len(chars)
    char_to_index = dict([(char, i) for (i, char) in enumerate(chars)])
    index_to_char = dict([(i, char) for (char, i) in char_to_index.items()])
    assert nchars == len(char_to_index) == len(index_to_char) > 1

    # define MCMC variable 
    priorvec = pymc.numpy.array([float(prefs[index_to_char[i]] * nchars * concentrationparam) for i in range(nchars)]) 
    initial_pi_incomplete = pymc.numpy.array([float(prefs[index_to_char[i]]) for i in range(nchars - 1)]) # incomplete Dirichlet 
    @pymc.stochastic(plot=False, name='prefs_posterior')
    def prefs_posterior(value=initial_pi_incomplete):
        """Posterior of preferences (value is incomplete Dirichlet vector)"""
        if any(value <= 0) or sum(value) >= 1:
            return -pymc.numpy.inf # outside Dirichlet support
        else:
            prior = pymc.distributions.dirichlet_like(value, priorvec)
            prefs_dict = dict([(index_to_char[i], pi_i) for (i, pi_i) in enumerate(value)])
            prefs_dict[index_to_char[nchars - 1]] = 1.0 - sum(value)
            tl.SetPreferences(prefs_dict, site)
            logl = tl.LogLikelihood()
            return prior + logl
    pi = pymc.CompletedDirichlet('pi', prefs_posterior)

    # run MCMC
    converged = False
    increasetry = 0
    while (not converged) and (increasetry <= nincreasetries):
        mcmc = pymc.MCMC([prefs_posterior, pi], verbose=max(0, verbose - 1))
        mcmc.use_step_method(DeltaExchange, prefs_posterior, verbose=max(0, verbose - 1))
        if verbose:
            print("\nBeginning MCMC of %d chains each of %d steps." % (nchains, nsteps))
        pi_means = {}
        for ichain in range(nchains):
            assert len(prefs_posterior.value) == nchars - 1
            prefs_posterior.value = pymc.numpy.random.dirichlet(priorvec)[ : -1] # initial random values
            mcmc.sample(iter=nsteps, burn=int(nsteps * burnfrac), progress_bar=(verbose >= 3))
            itrace = mcmc.trace('pi', chain=ichain)[:][:,0,:]
            assert itrace.shape[0] == nsteps - int(nsteps * burnfrac)
            pi_means[ichain] = sum(itrace) / float(itrace.shape[0])
            assert pymc.numpy.allclose(sum(pi_means[ichain]), 1.0)
        grstat = sum([r[0] for r in pymc.gelman_rubin(mcmc)['pi']]) / float(nchars)
        maxrmsdpi = max([math.sqrt(sum((pi_means[ichain] - pi_means[ichain + 1])**2)) for ichain in range(nchains - 1)])
        if (grstat < convergence) and (maxrmsdpi < (convergence - 1) / 2.0):
            converged = True
            mcmcstring = "MCMC converged (%d chains each with %d burn-in steps and %d sampling steps) with a Gelman-Rubin R of %.3f and a maximum RMS difference in preferences of %.3f." % (nchains, int(burnfrac * nsteps), nsteps - int(burnfrac * nsteps), grstat, maxrmsdpi)
            if verbose:
                print(mcmcstring)
        else:
            increasetry += 1
            mcmcstring = "MCMC FAILED to converge (%d chains each with %d burn-in steps and %d sampling steps) with a Gelman-Rubin R of %.3f and a maximum RMS difference in preferences of %.3f." % (nchains, int(burnfrac * nsteps), nsteps - int(burnfrac * nsteps), grstat, maxrmsdpi)
            if verbose:
                print(mcmcstring)
            if increasetry <= nincreasetries:
                nsteps = int(foldincrease * nsteps)
    traces = mcmc.trace('pi', chain=None)[:][:,0,:]
    assert traces.shape[0] == nchains * (nsteps - int(nsteps * burnfrac))
    pi_mean = sum(traces) / float(traces.shape[0])
    assert pymc.numpy.allclose(pi_mean, sum(pi_means.values()) / float(nchains))
    pi_mean = dict([(index_to_char[i], pi_i) for (i, pi_i) in enumerate(pi_mean)])
    return (pi_mean, mcmcstring, converged)


class DeltaExchange(pymc.Metropolis):
    """Implements delta exchange MCMC steps.

    Designed to operate on an incomplete Dirichlet distribution
    variable. Moves an amount from one element to another, effectively
    exchanging an amount delta between two elements of the distribution.
    Can propose moves that leave negative quantities.
    """
    def propose(self):
        """Moves random normal quantity from one element to another."""
        newvalue = self.stochastic.value.copy()
        n = len(newvalue)
        assert n >= 1
        i = pymc.numpy.random.randint(0, n + 1) # move from this element
        j = pymc.numpy.random.randint(0, n + 1) # move to this element
        while j == i:
            j = pymc.numpy.random.randint(0, n + 1)
        delta = pymc.numpy.random.normal(0, 0.05 * self.adaptive_scale_factor)
        if i < n: # not incompleted element
            newvalue[i] -= delta
        if j < n:
            newvalue[j] += delta
        self.stochastic.value = newvalue


if __name__ == '__main__':
    import doctest
    doctest.testmod()
