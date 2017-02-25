"""Tests branch length optimization by `TreeLikelihood`.

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


class test_BrLenOptimize_ExpCM(unittest.TestCase):
    """`TreeLikelihood` branch length optimization for `ExpCM`."""

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

        # amino-acid preferences
        self.nsites = 50
        prefs = []
        minpref = 0.02
        g = scipy.random.dirichlet([5] * N_NT)
        for r in range(self.nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))

        # simulate alignment with pyvolve
        pyvolvetree = pyvolve.read_tree(tree=self.newick)
        self.nseqs = self.tree.count_terminals()
        expcm = phydmslib.models.ExpCM(prefs)
        partitions = phydmslib.simulate.pyvolvePartitions(expcm)
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
        if self.MODEL == phydmslib.models.ExpCM:
            self.model = phydmslib.models.ExpCM(prefs)
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

    def test_Optimize(self):
        """Tests optimization of branch lengths."""
        tl = phydmslib.treelikelihood.TreeLikelihood(self.tree, 
                self.alignment, self.model, underflowfreq=self.underflowfreq)
        maxresult = tl.maximizeLikelihood(optimize_brlen=True)

        tl2 = phydmslib.treelikelihood.TreeLikelihood(self.tree, 
                self.alignment, self.model, underflowfreq=self.underflowfreq)
        maxresult = tl.maximizeLikelihood(optimize_brlen=False)

        self.assertTrue(tl.loglik > tl2.loglik)


class test_BrLenOptimize_ExpCM_gamma_omega(
        test_BrLenOptimize_ExpCM):
    MODEL = phydmslib.models.ExpCM
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedOmegaModel


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
