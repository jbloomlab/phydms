"""Tests `phydmslib.treelikelihood.TreeLikelihood`.

Written by Jesse Bloom.
"""

import os
import sys
import re
import unittest
import random
import copy
import scipy
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import Bio.Phylo
import phydmslib.treelikelihood
from phydmslib.constants import *


class test_TreeLikelihood(unittest.TestCase):
    """Tests `phydmslib.treelikelihood.TreeLikelihood` class."""

    def setUp(self):
        """Set up parameters for test."""
        random.seed(1)
        scipy.random.seed(1)

        # define tree and write image to a file
        self.newick = ('(($node_1=CAACAT$:0.1,$node_2=CAGCAG$:0.15)'
                       '$node_4=x$:0.15,$node_3=GAAAAG$:0.25)$node_5=y$:0.02;')
        tempfile = '_temp.tree'
        with open(tempfile, 'w') as f:
            f.write(self.newick)
        self.tree = Bio.Phylo.read(tempfile, 'newick')
        os.remove(tempfile)
        branch_labels = {} # branch annotations
        self.brlen = {}
        cladebyname = dict([(clade.name, clade) for clade in 
                self.tree.find_clades()])
        for (name, brlen) in re.findall(
                '(?P<name>\$node_\d\=[A-z]+\$)\:(?P<brlen>\d+\.\d+)', 
                self.newick):
            if name != self.tree.root.name:
                i = name.split('=')[0][-1] # node number
                branch_labels[cladebyname[name]] = "$t_{0}={1}$".format(
                        i, brlen)
                self.brlen[int(i)] = float(brlen)
        matplotlib.rc('text', usetex=True)
        Bio.Phylo.draw(self.tree, do_show=False, branch_labels=branch_labels)
        plt.axis('off')
        plt.savefig('test_treelikelihood_image.pdf')

        # define alignment
        self.nseqs = self.tree.count_terminals()
        self.nsites = 2
        self.alignment = []
        self.codons = {} # indexed by node, site, gives codon index
        for node in self.tree.get_terminals():
            seq = node.name.split('=')[1][ : -1]
            i = int(node.name.split('=')[0][-1]) # node number
            self.codons[i] = {}
            assert len(seq) == 3 * self.nsites
            self.alignment.append((node.name, seq))
            for r in range(self.nsites):
                codon = seq[3 * r : 3 * r + 3]
                self.codons[i][r] = CODON_TO_INDEX[codon]
        assert len(self.alignment) == self.nseqs

        # define model
        self.prefs = []
        minpref = 0.02
        for r in range(self.nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            self.prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
        phi = scipy.random.dirichlet([2] * N_NT)
        omega = 0.7
        kappa = 2.5
        beta = 1.6
        self.model = phydmslib.models.ExpCM(self.prefs, phi=phi, omega=omega,
                kappa=kappa, beta=beta)


    def test_InitializeTreeLikelihood(self):
        """Test that `TreeLikelihood` initializes properly."""
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree, self.alignment,
                self.model)
        self.assertTrue(tl.nsites == self.nsites)
        self.assertTrue(tl.nseqs == self.nseqs)
        self.assertTrue(tl.nnodes == len(tl.internalnodes) + self.nseqs)
        self.assertTrue(all([t > 0 for t in tl.t]))
        for n in tl.internalnodes:
            for descend in [tl.rdescend, tl.ldescend]:
                self.assertTrue(0 <= descend[n] < n, "{0}, {1}".format(n, descend[n]))

    def test_Likelihood(self):
        """Tests likelihood of `TreeLikelihood` object."""
        mus = [0.5, 1.5]
        partials_by_mu = {}
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
                M[node] = model.M(t)
            # compute partials at root node
            partials = scipy.zeros(shape=(self.nsites, N_CODON))
            for r in range(self.nsites):
                for y in range(N_CODON):
                    for x in range(N_CODON):
                        partials[r][y] += (M[3][r][y][self.codons[3][r]] *
                                M[4][r][y][x] * M[1][r][x][self.codons[1][r]]
                                * M[2][r][x][self.codons[2][r]])
            partials_by_mu[mu] = {'actual':tl.L[-1], 'expected':partials}

        for (i, mu1) in enumerate(mus):
            self.assertTrue(scipy.allclose(partials_by_mu[mu1]['actual'],
                    partials_by_mu[mu1]['expected']))
            for mu2 in mus[i + 1 : ]:
                self.assertFalse(scipy.allclose(partials_by_mu[mu1]['actual'],
                        partials_by_mu[mu2]['actual']))



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
