"""Tests correction for numerical underflow likelihood calculations.

Written by Jesse Bloom.
"""

import os
import sys
import scipy
import math
import unittest
import random
import phydmslib.models
import phydmslib.file_io
import phydmslib.treelikelihood
from phydmslib.constants import *
import Bio.Phylo


class test_underflow(unittest.TestCase):
    """Test underflow correction in `TreeLikelihood`."""

    def test_underflow(self):
        """Tests correction for numerical underflow."""

        scipy.random.seed(1)
        random.seed(1)

        treefile = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_tree.newick'))
        alignmentfile = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_alignment.fasta'))
        prefsfile = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_prefs.tsv'))

        alignment = phydmslib.file_io.ReadCodonAlignment(alignmentfile,
                checknewickvalid=True)
        tree = Bio.Phylo.read(treefile, 'newick')
        tree.root_at_midpoint()
        prefs = phydmslib.file_io.readPrefs(prefsfile, minpref=0.005)
        sites = sorted(prefs.keys())
        prefslist = [prefs[r] for r in sites]

        model = phydmslib.models.ExpCM(prefslist)
        distmodel = phydmslib.models.GammaDistributedOmegaModel(
                model, ncats=4)

        for (m, desc) in [(model, 'simple'), (distmodel, 'distribution')]:
            tl = phydmslib.treelikelihood.TreeLikelihood(tree,
                    alignment, m, underflowfreq=1)
            tl_nocorrection = phydmslib.treelikelihood.TreeLikelihood(
                    tree, alignment, m, underflowfreq=100000)

            self.assertTrue(scipy.allclose(tl.loglik, tl_nocorrection.loglik),
                    ("Log likelihoods differ with and without correction "
                    "for {0}: {1} versus {2}".format(desc, tl.loglik,
                    tl_nocorrection.loglik)))
            for (param, dl) in tl_nocorrection.dloglik.items():
                self.assertTrue(scipy.allclose(dl, tl.dloglik[param]),
                        ('dloglik differs with and without correction for {0}: '
                        '{1}, {2} versus {3}'.format(param, desc,
                        tl.dloglik[param], dl)))

            oldloglik = tl.loglik
            tl.dparamscurrent = False
            tl_nocorrection.dparamscurrent = False
            tl.dtcurrent = True
            tl_nocorrection.dtcurrent = True
            self.assertTrue(scipy.allclose(tl.loglik, tl_nocorrection.loglik),
                    ("Log likelihoods differ with and without correction "
                    "for {0}: {1} versus {2}".format(desc, tl.loglik,
                    tl_nocorrection.loglik)))
            self.assertTrue(scipy.allclose(tl.loglik, oldloglik))
            self.assertTrue(scipy.allclose(tl.dloglik_dt,
                    tl_nocorrection.dloglik_dt))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
