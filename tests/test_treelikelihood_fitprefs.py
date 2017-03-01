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
        kappa = 3.0
        omega = 3.0
        phi = scipy.random.dirichlet([5] * N_NT)
        self.model = phydmslib.models.ExpCM_fitprefs(self.prefs, 
                prior=None, kappa=kappa, omega=omega, phi=phi)
        self.realmodel = phydmslib.models.ExpCM(self.realprefs, 
                kappa=kappa, omega=omega, mu=10.0, phi=phi)

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

    def test_fitprefs_noprior(self):
        """Tests fitting of preferences with no prior."""
        firstloglik = None
        for seed in range(3):
            scipy.random.seed(seed)
            self.tl.paramsarray = scipy.random.uniform(0.1, 0.99, 
                    len(self.tl.paramsarray))
            if firstloglik:
                self.assertFalse(scipy.allclose(firstloglik, self.tl.loglik, 
                        atol=0.02), "loglik matches even with random pi")
            maxresult = self.tl.maximizeLikelihood()
            if firstloglik:
                self.assertTrue(scipy.allclose(firstloglik, self.tl.loglik,
                        atol=0.02), "loglik not the same for different starts")
            else:
                firstloglik = self.tl.loglik
            self.assertTrue(scipy.allclose(self.tl.model.origpi, self.model.pi))
            # ensure highest prefs are for amino acids with nonzero counts
            for r in range(self.tl.nsites):
                aas_with_counts = set([a for a in range(N_AA) if 
                        self.aacounts[r][a]])
                pref_sorted_aas = [tup[1] for tup in sorted(
                        [(self.tl.model.pi[r][a], a) for a in range(N_AA)],
                        reverse=True)]
                self.assertTrue(aas_with_counts == set(pref_sorted_aas[ : 
                        len(aas_with_counts)]), "top prefs not for amino "
                        "acids with counts.")


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
