"""Tests `phydmslib.treelikelihood.TreeLikelihood`.

Written by Jesse Bloom.
"""

import os
import sys
import re
import unittest
import random
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
        self.newick = ('(($node_1=CAA---$:0.1,$node_2=CAGCAG$:0.15)'
                       '$node_4=x$:0.15,$node_3=GAAAAG$:0.25)$node_5=y$:0.02;')
        tempfile = '_temp.tree'
        with open(tempfile, 'w') as f:
            f.write(self.newick)
        self.tree = Bio.Phylo.read(tempfile, 'newick')
        os.remove(tempfile)
        branch_labels = {} # branch annotations
        cladebyname = dict([(clade.name, clade) for clade in 
                self.tree.find_clades()])
        for (name, brlen) in re.findall(
                '(?P<name>\$node_\d\=[A-z]+\$)\:(?P<brlen>\d+\.\d+)', 
                self.newick):
            if name != self.tree.root.name:
                i = name.split('=')[0][-1] # node number
                branch_labels[cladebyname[name]] = "$t_{0}={1}$".format(
                        i, brlen)
        matplotlib.rc('text', usetex=True)
        Bio.Phylo.draw(self.tree, do_show=False, branch_labels=branch_labels)
        plt.axis('off')
        plt.savefig('test_treelikelihood_image.pdf')

        # define alignment
        self.nseqs = self.tree.count_terminals()
        self.nsites = 2
        self.alignment = []
        for node in self.tree.get_terminals():
            seq = node.name.split('=')[1][ : -1]
            assert len(seq) == 3 * self.nsites
            self.alignment.append((node.name, seq))
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


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
