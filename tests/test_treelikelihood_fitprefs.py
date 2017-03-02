"""Tests `TreeLikelihood` with `ExpCM_fitprefs`.

Written by Jesse Bloom."""


import random
import unittest
import copy
import os
import scipy
import scipy.linalg
import sympy
from phydmslib.constants import *
import phydmslib.models
import phydmslib.file_io
import phydmslib.simulate
import phydmslib.treelikelihood
import pyvolve


class test_TreeLikelihood_ExpCM_fitprefs(unittest.TestCase):
    """Test `ExpCM` with preferences as free parameters."""

    def setUp(self):
        """Set up for tests."""
        scipy.random.seed(1)
        random.seed(1)

        nsites = 1
        minpref = 0.001
        self.prefs = []
        self.realprefs = []
        for r in range(nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            self.prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
            scipy.random.shuffle(rprefs)
            self.realprefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
        self.kappa = 3.0
        self.omega = 3.0
        self.phi = scipy.random.dirichlet([5] * N_NT)
        self.model = phydmslib.models.ExpCM_fitprefs(self.prefs, 
                prior=None, kappa=self.kappa, omega=self.omega, phi=self.phi)
        self.realmodel = phydmslib.models.ExpCM(self.realprefs, 
                kappa=self.kappa, omega=self.omega, mu=10.0, phi=self.phi)

        treefile = 'NP_data/NP_tree.newick'
        self.tree = Bio.Phylo.read(treefile, 'newick')
        self.tree.root_at_midpoint()

        # simulate alignment using realmodel
        evolver = pyvolve.Evolver(
                partitions=phydmslib.simulate.pyvolvePartitions(self.realmodel),
                tree=pyvolve.read_tree(file=treefile))
        alignmentfile = '_temp_fitprefs_simulatedalignment.fasta'
        info = '_temp_info.txt'
        rates = '_temp_ratefile.txt'
        evolver(seqfile=alignmentfile, infofile=info, ratefile=rates)
        self.alignment = phydmslib.file_io.ReadCodonAlignment(alignmentfile, True)
        assert len(self.alignment[0][1]) == nsites * 3
        for f in [alignmentfile, info, rates]:
            os.remove(f)
        self.aacounts = dict([(r, dict([(a, 0) for a in range(N_AA)])) for
                r in range(nsites)])
        for (head, seq) in self.alignment:
            self.aacounts[r][CODON_TO_AA[CODON_TO_INDEX[seq]]] += 1

        self.tl = phydmslib.treelikelihood.TreeLikelihood(self.tree,
                self.alignment, self.model)

    def test_fitprefs_invquadratic_prior(self):
        """Tests fitting of preferences with invquadratic prior."""
        tls = {'no prior':copy.deepcopy(self.tl)}
        for (name, c1, c2) in [('weak prior', 100, 0.25),
                               ('strong C1 prior', 200, 0.25),
                               ('strong C2 prior', 100, 0.75)]:
            model = phydmslib.models.ExpCM_fitprefs(self.prefs,
                    prior=('invquadratic', c1, c2), kappa=self.kappa,
                    omega=self.omega, phi=self.phi)
            tls[name] = phydmslib.treelikelihood.TreeLikelihood(self.tree,
                    self.alignment, model)

        logliks = [tl.loglik for tl in tls.values()]
        self.assertTrue(all([scipy.allclose(logliks[0], logliki)
                for logliki in logliks]),
                "All loglik should be equal prior to optimization as pi "
                "starts at initial origpi values.")

        for (name, tl) in tls.items():
            maxresult = tl.maximizeLikelihood()

        self.assertTrue([tls['no prior'].loglik > tl.loglik for
                (name, tl) in tls.items() if name != 'no prior'],
                "Unconstrained model does not have the highest loglik.")

        self.assertTrue([tls['weak prior'].loglik > tl.loglik for
                (name, tl) in tls.items() if 'strong' in name],
                "Weak prior does not have higher loglik than strong prior.")

        pidiff2 = dict([(name, ((tl.model.pi - tl.model.origpi)**2).sum()) for
                (name, tl) in tls.items()])
        self.assertTrue([pidiff2['no prior'] > x for (name, x) in
                pidiff2.items() if name != 'no prior'])
        self.assertTrue([pidiff2['weak prior'] > x for (name, x) in
                pidiff2.items() if 'strong' in name])


    def test_fitprefs_noprior(self):
        """Tests fitting of preferences with no prior."""
        tl = copy.deepcopy(self.tl)
        firstloglik = None
        for seed in range(3):
            scipy.random.seed(seed)
            tl.paramsarray = scipy.random.uniform(0.1, 0.99, 
                    len(tl.paramsarray))
            if firstloglik:
                self.assertFalse(scipy.allclose(firstloglik, tl.loglik, 
                        atol=0.02), "loglik matches even with random pi")
            maxresult = tl.maximizeLikelihood()
            if firstloglik:
                self.assertTrue(scipy.allclose(firstloglik, tl.loglik,
                        atol=0.02), "loglik not the same for different starts")
            else:
                firstloglik = tl.loglik
            self.assertTrue(scipy.allclose(tl.model.origpi, self.model.pi))
            # ensure highest prefs are for amino acids with nonzero counts
            for r in range(tl.nsites):
                aas_with_counts = set([a for a in range(N_AA) if 
                        self.aacounts[r][a]])
                pref_sorted_aas = [tup[1] for tup in sorted(
                        [(tl.model.pi[r][a], a) for a in range(N_AA)],
                        reverse=True)]
                self.assertTrue(aas_with_counts <= set(pref_sorted_aas[ : 
                        len(aas_with_counts) + 1]), "top prefs not for amino "
                        "acids with counts:\naas_with_counts: {0}\n"
                        "pref_sorted_aas: {1}".format(aas_with_counts,
                        pref_sorted_aas))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
