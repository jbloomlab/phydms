"""Module for MCMC; uses ``pymc`` (tested with version 2.3.4)."""


import math
import random
import pymc

def AdjustToMinValue(vec, minvalue):
    """Adjusts incomplete Dirichlet vector *vec* so all elements > *minvalue*"""
    nadjusts = 0
    vec = pymc.numpy.ravel(vec)
    while any(vec <= minvalue) or sum(vec) >= 1 - minvalue:
        maxindex = pymc.numpy.argmax(vec)
        if 1 - sum(vec) > vec[maxindex]:
            maxindex = -1
        else:
            vec[maxindex] -= minvalue
        minindex = pymc.numpy.argmin(vec)
        if 1 - sum(vec) < vec[minindex]:
            minindex = -1
        else:
            vec[minindex] += minvalue
        assert maxindex != minindex
        nadjusts += 1
        assert nadjusts < 2 * len(vec), "problem adjusting"
    return vec


def PrefsMCMC(tl, prefs, site, concentrationparam, seed, verbose=0, nchains=2, nsteps=1e4, burnfrac=0.2, nincreasetries=4, convergence=1.1, minvalue=1.0e-3, use_delta_exchange=False):
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

    *nchains*, *burnfrac*, *nsteps*, *nincreasetries*, and 
    *convergence* specify how the MCMC is performed. The MCMC
    is considered to converge if the chains (*nchains* must be >= 2)
    have a Gelman-Rubin R that is less than *convergence* **and**
    if the max root-mean-square difference in the inferred preferences
    between chains is less than *(convergence - 1.0) / 2*. Each chain begins
    with *nsteps X burnfrac* steps. Then we perform a chain of *nsteps*
    steps. If convergence fails, we add *nsteps* additional steps to
    the chain *nincreasetries* more times. Each chain is seeded with
    preferences drawn randomly from the prior.

    We do not allow any element of the preference vector to become
    smaller than *minvalue*.

    *use_delta_exchange* specifies that we use the delta exchange
    MCMC operator rather than the ``pymc`` default.

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
    assert nincreasetries >= 0
    assert concentrationparam > 0, "concentratioparam must be > 0"
    assert all([pi > 0 for pi in prefs.values()]), "The preference must be > 0 for all sites"
    assert abs(sum(prefs.values()) - 1.0) < 1.0e-5, "The sum of the preferences must be one"
    assert isinstance(site, int) and 1 <= site <= tl.NSites(), "site of %d is not in the tree likelihood object" % site
    assert set(prefs.keys()) == set(tl.GetPreferences(site).keys()), "prefs does not have keys for the same characters as the tree likelihood object"

    # seed random number generator
    pymc.numpy.random.seed(seed)

    # dicts to convert from prefs chars (e.g. amino acids) to vector indices
    # order from smallest to largest; this makes last element of vector the
    # incompleted item in Dirichlet which can help with convergence; when steps
    # are specified for other variables last one will have biggest change
    decorated_chars = [(ipi, ichar) for (ichar, ipi) in prefs.items()]
    decorated_chars.sort()
    chars = [tup[1] for tup in decorated_chars]
    nchars = len(chars)
    char_to_index = dict([(char, i) for (i, char) in enumerate(chars)])
    index_to_char = dict([(i, char) for (char, i) in char_to_index.items()])
    assert nchars == len(char_to_index) == len(index_to_char) > 1

    # define MCMC variables 
    priorvec = pymc.numpy.array([float(prefs[index_to_char[i]] * nchars * concentrationparam) for i in range(nchars)]) 
    initial_pi_incomplete = AdjustToMinValue(pymc.numpy.array([float(prefs[index_to_char[i]]) for i in range(nchars - 1)]), minvalue) # incomplete Dirichlet 

    @pymc.stochastic(plot=False, name='prefs_posterior')
    def prefs_posterior(value=initial_pi_incomplete):
        """Posterior of preferences (value is incomplete Dirichlet vector)"""
        if any(value <= minvalue) or sum(value) >= 1 - minvalue:
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
    tune_interval = 250
    converged = False
    increasetry = 0
    traces = {} # keyed by chain
    mcmc = {} # keyed by chain
    while (not converged) and (increasetry <= nincreasetries):
        increasetry += 1
        for ichain in range(1, nchains + 1):
            assert len(prefs_posterior.value) == nchars - 1
            if ichain not in mcmc:
                mcmc[ichain] = pymc.MCMC([prefs_posterior, pi], verbose=max(0, verbose - 1))
                if use_delta_exchange:
                    mcmc[ichain].use_step_method(DeltaExchange, prefs_posterior, verbose=max(0, verbose - 1))
                else:
                    mcmc[ichain].use_step_method(pymc.Metropolis, prefs_posterior, verbose=max(0, verbose - 1))
                prefs_posterior.value = AdjustToMinValue(pymc.numpy.random.dirichlet(priorvec)[ : -1], 0.01) # initial value drawn from posterior to start, don't let any element be too small (hence 0.01 adjustment)
                if verbose:
                    print("Beginning burn-in of %d steps for chain %d, starting with %s" % (int(nsteps * burnfrac), ichain, str(prefs_posterior.value)))
                mcmc[ichain].sample(iter=int(nsteps * burnfrac), burn=int(nsteps * burnfrac), progress_bar=(verbose >= 3), tune_interval=tune_interval) # burn in
                assert len(mcmc[ichain].trace('pi', chain=None)[:][:,0,:]) == 1, "Found traces after burnin: %s" % str(mcmc[ichain].trace('pi', chain=None)[:][:,0,:])
                if verbose:
                    print("Completed burn-in for chain %d, values now %s" % (ichain, str(prefs_posterior.value)))
                mcmc[ichain].save_state()
            mcmc[ichain].restore_sampler_state()
            if verbose:
                print("Now sampling %d steps from chain %d for the %d time, starting with %s" % (nsteps, ichain, increasetry, str(prefs_posterior.value)))
            mcmc[ichain].sample(iter=nsteps, burn=0, progress_bar=(verbose >= 3), tune_interval=tune_interval)
            mcmc[ichain].save_state()
            traces[ichain] = mcmc[ichain].trace('pi', chain=None)[:][:,0,:][1 : ] # first sample is burn-in empty vector, and so removed
            assert len(traces[ichain]) == increasetry * nsteps, "Failed to find expected chain length after %d increase samplings: %s, length %d" % (increasetry, str(traces[ichain]), len(traces[ichain]))
            if verbose:
                print("Finished sampling from chain %d for the %d time, values now %s" % (ichain, increasetry, str(prefs_posterior.value)))
        pi_means = dict([(ichain, sum(itrace) / float(itrace.shape[0])) for (ichain, itrace) in traces.items()])
        assert all([pymc.numpy.allclose(sum(imean), 1.0) for imean in pi_means.values()]), str(pi_means)
        maxrmsdpi = max([math.sqrt(sum((pi_means[ichain] - pi_means[ichain + 1])**2)) for ichain in range(1, nchains)])
        grstatvec = [r for r in pymc.gelman_rubin(pymc.numpy.array([itrace for itrace in traces.values()]))]
        grstat = sum([r for r in grstatvec]) / float(nchars)
        if verbose:
            print("Examined convergence. The max RMS diff pref is %g, the average Gelman-Rubin statistic is %g, and the Gelman-Rubin vector is %s" % (maxrmsdpi, grstat, grstatvec))
        if (grstat < convergence) and (maxrmsdpi < (convergence - 1) / 2.0):
            converged = True
        else:
            converged = False
        mcmcstring = "MCMC %s (%d chains each with %d sampling steps) with a Gelman-Rubin R of %.3f and a maximum RMS difference in preferences of %.3f." % ({True:'converged', False:'FAILED to converge'}[converged], nchains, nsteps * increasetry, grstat, maxrmsdpi)
        if verbose:
            print(mcmcstring)
    assert len(pi_means) == nchains
    pi_mean = sum(pi_means.values()) / float(len(pi_means))
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
