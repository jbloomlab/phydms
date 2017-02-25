"""Tests branch length derivatives by `TreeLikelihood`.

Written by Jesse Bloom.
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


class test_BrLenDerivatives_ExpCM(unittest.TestCase):
    """`TreeLikelihood` branch length derivatives for `ExpCM`."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM
    DISTRIBUTIONMODEL = None

    def setUp(self):
        """Set up parameters for test."""
        random.seed(1)
        scipy.random.seed(1)

        self.underflowfreq = 1

        # define tree 
        self.newick = ('((node1:0.2,node2:0.3)node4:0.3,node3:0.5)node5:0.04;')
        tempfile = '_temp.tree'
        with open(tempfile, 'w') as f:
            f.write(self.newick)
        self.tree = Bio.Phylo.read(tempfile, 'newick')
        os.remove(tempfile)

        # simulate alignment with pyvolve
        pyvolvetree = pyvolve.read_tree(tree=self.newick)
        self.nsites = 50
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
        else:
            raise ValueError("Invalid DISTRIBUTIONMODEL: {0}".format(
                    self.DISTRIBUTIONMODEL))

    def test_Initialize(self):
        """Test that `TreeLikelihood` initializes properly."""
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree, 
                self.alignment, self.model, underflowfreq=self.underflowfreq,
                dparamscurrent=False, dtcurrent=True)
        self.assertEqual(tl.dloglik_dt.shape, tl.t.shape)

    def test_dM_dt(self):
        """Tests model `dM` with respect to `t`."""
        scipy.random.seed(1)
        random.seed(1)
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree, 
                self.alignment, self.model, underflowfreq=self.underflowfreq,
                dparamscurrent=False, dtcurrent=True)

        def func(t, k, r, x, y):
            return tl._M(k, t[0], None)[r][x][y]

        def dfunc(t, k, r, x, y):
            return tl._dM(k, t[0], 't', None)[r][x][y]

        for r in range(self.nsites):
            k = random.choice(tl._catindices)
            x = random.randint(0, N_CODON - 1)
            y = random.randint(0, N_CODON - 1)
            t = scipy.array([random.uniform(0.05, 0.5)])
            diff = scipy.optimize.check_grad(func, dfunc, t, k, r, x, y)
            self.assertTrue(abs(diff) < 1e-4, "diff = {0}".format(diff))

    def test_AdjustBrLen(self):
        """Tests adjusting branch lengths."""
        scipy.random.seed(1)
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree, 
                self.alignment, self.model, underflowfreq=self.underflowfreq,
                dparamscurrent=False, dtcurrent=True)
        loglik1 = tl.loglik
        tl.t = 2 * tl.t
        self.assertFalse(scipy.allclose(loglik1, tl.loglik))
        tl.t = tl.t / 2
        self.assertTrue(scipy.allclose(loglik1, tl.loglik))
        tl.t = tl.t * scipy.random.uniform(0.1, 1.5, tl.t.shape)
        self.assertFalse(scipy.allclose(loglik1, tl.loglik))
        loglik2 = tl.loglik
        t = tl.t
        t[1] *= 3
        tl.t = t
        self.assertFalse(scipy.allclose(loglik2, tl.loglik))
        t = tl.t
        t[1] /= 3
        tl.t = t
        self.assertTrue(scipy.allclose(loglik2, tl.loglik))

    def test_BrLenDerivatives(self):
        """Tests derivatives of branch lengths."""
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree, 
                self.alignment, self.model, underflowfreq=self.underflowfreq,
                dparamscurrent=False, dtcurrent=True)

        def func(t, n):
            tx = tl.t
            tx[n] = t[0]
            tl.t = tx
            return tl.loglik

        def dfunc(t, n):
            tx = tl.t
            tx[n] = t[0]
            tl.t = tx
            return tl.dloglik_dt[n]

        for n in range(len(tl.t)):
            diff = scipy.optimize.check_grad(func, dfunc, 
                    scipy.array([tl.t[n]]), n)
            self.assertTrue(diff < 1e-5, diff)


    def test_dtcurrent(self):
        """Tests use of `dtcurrent` attribute."""
        tl1 = phydmslib.treelikelihood.TreeLikelihood(self.tree, 
                self.alignment, self.model, underflowfreq=self.underflowfreq,
                dparamscurrent=False, dtcurrent=True)

        tl2 = phydmslib.treelikelihood.TreeLikelihood(self.tree, 
                self.alignment, self.model, underflowfreq=self.underflowfreq,
                dparamscurrent=True, dtcurrent=False)

        self.assertTrue(scipy.allclose(tl1.loglik, tl2.loglik))

        with self.assertRaises(Exception) as context:
            tl2.dtcurrent = True

        with self.assertRaises(Exception) as context:
            tl1.dparamscurrent = True

        with self.assertRaises(Exception) as context:
            x = tl2.dloglik_dt

        with self.assertRaises(Exception) as context:
            x = tl1.dloglik

        tl1.t = 1.1 * tl1.t
        tl2.t = 1.1 * tl2.t
        tl2.dparamscurrent = False
        tl2.dtcurrent = True
        self.assertTrue(scipy.allclose(tl1.dloglik_dt, tl2.dloglik_dt))
        self.assertTrue(scipy.allclose(tl1.loglik, tl2.loglik))

        tl1.updateParams({'kappa':3.15})
        tl2.updateParams({'kappa':3.15})
        tl2.dtcurrent = False
        tl2.dparamscurrent = True
        tl1.dtcurrent = False
        tl1.dparamscurrent = True
        for (param, deriv) in tl1.dloglik.items():
            self.assertTrue(scipy.allclose(deriv, tl2.dloglik[param]))
        self.assertTrue(scipy.allclose(tl1.loglik, tl2.loglik))


class test_BrLenDerivatives_ExpCM_empirical_phi(test_BrLenDerivatives_ExpCM):
    MODEL = phydmslib.models.ExpCM_empirical_phi


class test_BrLenDerivativs_ExpCM_empirical_phi_divpressure(test_BrLenDerivatives_ExpCM):
    MODEL = phydmslib.models.ExpCM_empirical_phi_divpressure


class test_BrLenDerivatives_YNGKP_M0(test_BrLenDerivatives_ExpCM):
    MODEL = phydmslib.models.YNGKP_M0


class test_BrLenDerivatives_YNGKP_M5(test_BrLenDerivatives_ExpCM):
    MODEL = phydmslib.models.YNGKP_M0
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedOmegaModel


class test_BrLenDerivatives_ExpCM_gamma_omega(
        test_BrLenDerivatives_ExpCM):
    MODEL = phydmslib.models.ExpCM
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedOmegaModel


class test_BrLenDerivatives_ExpCM_empirical_phi_gamma_omega(
        test_BrLenDerivatives_ExpCM):
    MODEL = phydmslib.models.ExpCM_empirical_phi
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedOmegaModel


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
