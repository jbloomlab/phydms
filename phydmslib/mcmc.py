"""Module for MCMC."""


import math
import math
import random
import pymc


def PrefsMCMC(tl, prefs, site, concentrationparam):
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

    The return value is *(meanprefs, mcmcstring)* where:

        - *meanprefs* is a dictionary with the same keys as *prefs*,
          and with values giving the posterior mean preference for
          that character.

        - *mcmcstring* is a string that summarizes some information
          about the MCMC and its convergence.
    """
    assert concentrationparam > 0, "concentratioparam must be > 0"
    assert all([pi > 0 for pi in prefs.values()]), "The preference must be > 0 for all sites"
    assert abs(sum(pi.values()) - 1.0) < 1.0e-5, "The sum of the preferences must be one"
    assert isinstance(site, int) and 1 <= site <= tl.NSites(), "site of %d is not in the tree likelihood object" % site
    assert set(prefs.keys()) == set(tl.GetPreferences().keys()), "prefs does not have keys for the same characters as the tree likelihood object"


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
        i = random.randint(0, n) # move from this element
        j = random.randint(0, n) # move to this element
        while j == i:
            j = random.randint(0, n)
        delta = random.normalvariate(0, 0.1 * self.adaptive_scale_factor)
        if i < n: # not incompleted element
            newvalue[i] -= delta
        if j < n:
            newvalue[j] += delta
        self.stochastic.value = newvalue
        #print "For %s, moved %g from %d to %d (adaptive_scale_factor = %g)" % (self.stochastic.__name__, delta, i, j, self.adaptive_scale_factor)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
