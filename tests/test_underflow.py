"""Tests correction for numerical underflow in `TreeLikelihood`.

Written by Jesse Bloom.
"""

import os
import sys
import scipy
import math
import unittest
import random
import io
import copy
import phydmslib.models
import phydmslib.file_io
import phydmslib.treelikelihood
from phydmslib.constants import *
import Bio.SeqIO
import Bio.Phylo
import pyvolve



class test_underflow(unittest.TestCase):
    """Tests correction for numerical underflow."""

    def test_underflow(self):
        """Tests correction for numerical underflow."""

        scipy.random.seed(1)
        random.seed(1)

        alignmentfile = 'NP_data/NP_alignment.fasta'
        treefile = 'NP_data/NP_tree.newick'
        prefsfile = 'NP_data/NP_prefs.csv'

        alignment = phydmslib.file_io.ReadCodonAlignment(alignmentfile,
                checknewickvalid=True)
        tree = Bio.Phylo.read(treefile, 'newick')
        tree.root_at_midpoint()
        prefs = phydmslib.file_io.readPrefs(prefsfile, minpref=0.005)
        sites = sorted(prefs.keys())
        prefslist = [prefs[r] for r in sites]

        model = phydmslib.models.ExpCM(prefslist)

        tl = phydmslib.treelikelihood.TreeLikelihood(tree, alignment, model,
                underflowfreq=5)
        tl_nocorrection = phydmslib.treelikelihood.TreeLikelihood(tree, 
                alignment, model, underflowfreq=10000)

        self.assertTrue(scipy.allclose(tl.loglik, tl_nocorrection.loglik),
                ("Log likelihoods differ with and without correction: "
                "{0} versus {1}".format(tl.loglik, tl_nocorrection.loglik)))
        for (param, dl) in tl_nocorrection.dloglik.items():
            self.assertTrue(scipy.allclose(dl, tl.dloglik[param]),
                    ('dloglik differs with and without correction for {0}: '
                    '{1} versus {2}'.format(param, tl.dloglik[param], dl)))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
