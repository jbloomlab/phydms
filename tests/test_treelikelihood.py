"""Tests `TreeLikelihood` for various models.

Written by Jesse Bloom and Sarah Hilton.
"""

import os
import sys
import re
import math
import unittest
import random
import copy
import scipy
import scipy.optimize
import Bio.Phylo
import phydmslib.models
import phydmslib.treelikelihood
import phydmslib.simulate
from phydmslib.constants import *
import pyvolve


class test_TreeLikelihood_ExpCM(unittest.TestCase):
    """Tests `TreeLikelihood` for `ExpCM` model."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM
    DISTRIBUTIONMODEL = None

    def setUp(self):
        """Set up parameters for test."""
        random.seed(1)
        scipy.random.seed(1)

        # define tree
        self.newick = ('((node1:0.2,node2:0.3)node4:0.3,node3:0.5)node5:0.04;')
        tempfile = '_temp.tree'
        with open(tempfile, 'w') as f:
            f.write(self.newick)
        self.tree = Bio.Phylo.read(tempfile, 'newick')
        os.remove(tempfile)
        self.brlen = {}
        for (name, brlen) in re.findall(
                '(?P<name>node\d):(?P<brlen>\d+\.\d+)', self.newick):
            if name != self.tree.root.name:
                i = name[-1] # node number
                self.brlen[int(i)] = float(brlen)

        # simulate alignment with pyvolve
        pyvolvetree = pyvolve.read_tree(tree=self.newick)
        self.nsites = 60
        self.nseqs = self.tree.count_terminals()
        e_pw = scipy.ndarray((3, N_NT), dtype='float')
        e_pw.fill(0.25)
        yngkp_m0 = phydmslib.models.YNGKP_M0(e_pw, self.nsites)
        partitions = phydmslib.simulate.pyvolvePartitions(yngkp_m0)
        alignment = '_temp_simulatedalignment.fasta'
        info = '_temp_info.txt'
        rates = '_temp_ratefile.txt'
        evolver = pyvolve.Evolver(partitions=partitions, tree=pyvolvetree)
        evolver(seqfile=alignment, infofile=info, ratefile=rates)
        self.alignment = [(s.description, str(s.seq)) for s in Bio.SeqIO.parse(
                alignment, 'fasta')]
        for f in [alignment, info, rates]:
            os.remove(f)
        assert len(self.alignment[0][1]) == self.nsites * 3
        assert len(self.alignment) == self.nseqs
        self.codons = {} # indexed by node, site, gives codon index
        for node in self.tree.get_terminals():
            node = node.name
            i = int(node[-1])
            self.codons[i] = {}
            seq = [seq for (head, seq) in self.alignment if node == head][0]
            for r in range(self.nsites):
                codon = seq[3 * r : 3 * r + 3]
                self.codons[i][r] = CODON_TO_INDEX[codon]

        # define model
        prefs = []
        minpref = 0.02
        g = scipy.random.dirichlet([5] * N_NT)
        for r in range(self.nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
        if self.MODEL == phydmslib.models.ExpCM:
            self.model = phydmslib.models.ExpCM(prefs)
        elif self.MODEL == phydmslib.models.ExpCM_empirical_phi:
            self.model = phydmslib.models.ExpCM_empirical_phi(prefs, g)
        elif self.MODEL == phydmslib.models.ExpCM_empirical_phi_divpressure:
            divpressure = scipy.random.uniform(-1, 5, self.nsites)
            divpressure /= max(abs(divpressure))
            self.model = phydmslib.models.ExpCM_empirical_phi_divpressure(
                    prefs, g, divpressure)
        elif self.MODEL == phydmslib.models.YNGKP_M0:
            e_pw = scipy.random.uniform(0.2, 0.8, size=(3, N_NT))
            e_pw = e_pw / e_pw.sum(axis=1, keepdims=True)
            self.model = phydmslib.models.YNGKP_M0(e_pw, self.nsites)
        else:
            raise ValueError("Invalid MODEL: {0}".format(self.MODEL))

        if self.DISTRIBUTIONMODEL is None:
            pass
        elif (self.DISTRIBUTIONMODEL ==
                phydmslib.models.GammaDistributedOmegaModel):
            self.model = self.DISTRIBUTIONMODEL(self.model, ncats=4)
        elif (self.DISTRIBUTIONMODEL ==
                phydmslib.models.GammaDistributedBetaModel):
            self.model = self.DISTRIBUTIONMODEL(self.model, ncats=4)
        else:
            raise ValueError("Invalid DISTRIBUTIONMODEL: {0}".format(
                    self.DISTRIBUTIONMODEL))

    def test_Initialize(self):
        """Test that initializes properly."""
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree,
                self.alignment, self.model)
        self.assertTrue(tl.nsites == self.nsites)
        self.assertTrue(tl.nseqs == self.nseqs)
        self.assertTrue(tl.nnodes == tl.ninternal + tl.ntips)
        self.assertTrue(tl.ntips == self.nseqs)
        self.assertTrue(all([t > 0 for t in tl.t]))
        for n in range(tl.ntips, tl.nnodes):
            for descend in [tl.rdescend, tl.ldescend]:
                i = n - tl.ntips
                self.assertTrue(0 <= descend[i] < n,
                        "{0}, {1}".format(n, descend[i]))
        self.assertTrue(tl.nsites == len(tl.siteloglik))

    def test_paramsarray(self):
        """Tests params array setting and getting."""
        modelparams = self.getModelParams(seed=1)
        model = copy.deepcopy(self.model)
        model.updateParams(modelparams)
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree,
                self.alignment, model)
        logl = tl.loglik
        paramsarray = tl.paramsarray
        nparams = len(paramsarray)
        self.assertTrue(nparams == sum(map(lambda x: (1 if isinstance(x, float)
                else len(x)), modelparams.values())))
        # set to new value, make sure have changed
        tl.paramsarray = scipy.array([random.uniform(0.7, 0.8) for i in
                range(nparams)])
        for (param, value) in modelparams.items():
            self.assertFalse(scipy.allclose(value, getattr(tl.model, param)))
        self.assertFalse(scipy.allclose(logl, tl.loglik))
        # re-set to old value, make sure return to original values
        tl.paramsarray = copy.deepcopy(paramsarray)
        self.assertTrue(scipy.allclose(logl, tl.loglik))
        for (param, value) in modelparams.items():
            self.assertTrue(scipy.allclose(value, getattr(tl.model, param)))

    def test_Likelihood(self):
        """Tests likelihood."""
        if self.DISTRIBUTIONMODEL:
            return # test doesn't work for DistributionModel
        mus = [0.5, 1.5]
        partials_by_mu = {}
        siteloglik_by_mu = {}
        loglik_by_mu = {}
        for mu in mus:
            model = copy.deepcopy(self.model)
            model.updateParams({'mu':mu})
            tl = phydmslib.treelikelihood.TreeLikelihood(self.tree,
                    self.alignment, model)
            # Here we are doing the multiplication hand-coded for the
            # tree defined in `setUp`. This calculation would be wrong
            # if the tree in `setUp` were to be changed.
            M = {}
            for (node, t) in self.brlen.items():
                M[node] = model.M(t / model.branchScale)
            # compute partials at root node
            partials = scipy.zeros(shape=(self.nsites, N_CODON))
            siteloglik = scipy.zeros(shape=(self.nsites,))
            loglik = 0.0
            for r in range(self.nsites):
                for y in range(N_CODON):
                    for x in range(N_CODON):
                        partials[r][y] += (M[3][r][y][self.codons[3][r]] *
                                M[4][r][y][x] * M[1][r][x][self.codons[1][r]]
                                * M[2][r][x][self.codons[2][r]])
                    siteloglik[r] += partials[r][y] * model.stationarystate[r, y]
                siteloglik[r] = math.log(siteloglik[r])
                loglik += siteloglik[r]
            rootnode = tl.nnodes - 1
            partials_by_mu[mu] = {'actual':tl.L[rootnode], 'expected':partials}
            siteloglik_by_mu[mu] = {'actual':tl.siteloglik, 'expected':siteloglik}
            loglik_by_mu[mu] = {'actual':tl.loglik, 'expected':loglik}

        for (i, mu1) in enumerate(mus):
            for (name, d) in [('partials', partials_by_mu),
                                      ('siteloglik', siteloglik_by_mu),
                                      ('loglik', loglik_by_mu)]:
                self.assertTrue(scipy.allclose(d[mu1]['actual'],
                        d[mu1]['expected']), "Mismatch: {0}".format(name))

    def test_LikelihoodDerivativesModelParams(self):
        """Test derivatives of with respect to model params."""
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree,
                self.alignment, self.model)

        for itest in range(2):
            modelparams = self.getModelParams(seed=itest)
            tl.updateParams(modelparams)

            def func(x, i):
                y = tl.paramsarray
                y[i] = x[0]
                tl.paramsarray = y
                return tl.loglik
            def dfunc(x, i):
                y = tl.paramsarray
                y[i] = x[0]
                tl.paramsarray = y
                return tl.dloglikarray[i]
            for iparam in range(len(tl.paramsarray)):
                diff = scipy.optimize.check_grad(func, dfunc,
                        scipy.array([tl.paramsarray[iparam]]), iparam)
                self.assertTrue(diff < 2e-3, ("{0}: diff {1}, value "
                        "{2}, deriv {3}").format(tl._index_to_param[iparam],
                        diff, tl.paramsarray[iparam],
                        dfunc([tl.paramsarray[iparam]], iparam)))

    def test_MaximizeLikelihood(self):
        """Tests maximization likelihood.

        Make sure it gives the same value for several starting points."""
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree,
                self.alignment, self.model)

        logliks = []
        paramsarrays = []
        for itest in range(3):
            modelparams = self.getModelParams(itest)
            tl.updateParams(modelparams)
            startloglik = tl.loglik
            result = tl.maximizeLikelihood()
            self.assertTrue(tl.loglik > startloglik, "no loglik increase: "
                    "start = {0}, end = {1}".format(startloglik, tl.loglik))
            for (otherloglik, otherparams) in zip(logliks, paramsarrays):
                self.assertTrue(scipy.allclose(tl.loglik, otherloglik,
                        atol=1e-3, rtol=1e-3),
                        "Large difference in loglik: {0} vs {1}".format(
                        otherloglik, tl.loglik))
                self.assertTrue(scipy.allclose(tl.paramsarray, otherparams,
                        atol=1e-2, rtol=1e-1),
                        "Large difference in paramsarray: {0} vs {1}, {2}".format(
                        otherparams, tl.paramsarray, self.model))
            logliks.append(tl.loglik)
            paramsarrays.append(tl.paramsarray)

    def getModelParams(self, seed):
        """Dict of random model params for random number seed `seed`."""
        random.seed(seed)
        scipy.random.seed(seed)
        modelparams = {}
        for param in self.model.freeparams:
            pvalue = getattr(self.model, param)
            if isinstance(pvalue, float):
                modelparams[param] = random.uniform(0.5, 1.5)
            elif isinstance(pvalue, scipy.ndarray) and pvalue.ndim == 1:
                modelparams[param] = scipy.random.dirichlet([6] * len(pvalue))
            else:
                raise ValueError("Can't handle {0}".format(param))
        return modelparams


class test_TreeLikelihood_ExpCM_empirical_phi(test_TreeLikelihood_ExpCM):
    """Tests `TreeLikelihood` for `ExpCM_empirical_phi`."""
    MODEL = phydmslib.models.ExpCM_empirical_phi


class test_TreeLikelihood_ExpCM_empirical_phi_divpressure(test_TreeLikelihood_ExpCM):
    """Tests `TreeLikelihood` for `ExpCM_empirical_phi_divpressure`."""
    MODEL = phydmslib.models.ExpCM_empirical_phi_divpressure


class test_TreeLikelihood_YNGKP_M0(test_TreeLikelihood_ExpCM):
    """Tests `TreeLikelihood` for `YNGKP_M0`."""
    MODEL = phydmslib.models.YNGKP_M0


class test_TreeLikelihood_YNGKP_M5(test_TreeLikelihood_ExpCM):
    """Tests for `YNGKP_M5`."""
    MODEL = phydmslib.models.YNGKP_M0
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedOmegaModel

class test_TreeLikelihood_ExpCM_empirical_phi_gamma_omega(
        test_TreeLikelihood_ExpCM):
    """Tests for `ExpCM_empirical_phi` with gamma omega."""
    MODEL = phydmslib.models.ExpCM_empirical_phi
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedOmegaModel

class test_TreeLikelihood_ExpCM_empirical_phi_gamma_beta(
        test_TreeLikelihood_ExpCM):
    """Tests for `ExpCM_empirical_phi` with gamma beta."""
    MODEL = phydmslib.models.ExpCM_empirical_phi
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedBetaModel

class test_TreeLikelihood_ExpCM_gamma_beta(
        test_TreeLikelihood_ExpCM):
    """Tests for `ExpCM` with gamma beta."""
    MODEL = phydmslib.models.ExpCM
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedBetaModel


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
