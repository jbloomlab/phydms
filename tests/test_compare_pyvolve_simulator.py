"""Tests whether or not `pyvolve` and `Simulator` give similar results.

Makes sure `pyvolve` and `Simulator` give similar distributions of
alignment summary statistics.

Written by Sarah Hilton and Jesse Bloom.
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
import phydmslib.treelikelihood
import phydmslib.simulate
from phydmslib.constants import *
from phydmslib.file_io import ReadCodonAlignment
import Bio.SeqIO
import Bio.Phylo
import pyvolve
import itertools
from scipy.special import comb
from scipy.stats import entropy
import pandas as pd


class test_compare_ExpCM(unittest.TestCase):
    """Tests `pyvolve` and `Simulator` on `ExpCM`"""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi

    def setUp(self):
        """Set up model, tree, and simulate alignments."""
        scipy.random.seed(1)
        # number of simulation replicates
        self.nsim = 100
        # define model
        self.nsites = 10
        prefs = []
        minpref = 0.01
        for r in range(self.nsites):
            rprefs = scipy.random.dirichlet([1] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
        kappa = 4.2
        omega = 0.4
        beta = 1.5
        mu = 0.3
        if self.MODEL == phydmslib.models.ExpCM:
            phi = scipy.random.dirichlet([7] * N_NT)
            self.model = phydmslib.models.ExpCM(prefs, kappa=kappa, omega=omega,
                                                beta=beta, mu=mu, phi=phi,
                                                freeparams=['mu'])
        elif self.MODEL == phydmslib.models.ExpCM_empirical_phi:
            g = scipy.random.dirichlet([7] * N_NT)
            self.model = phydmslib.models.ExpCM_empirical_phi(prefs, g,
                                                              kappa=kappa,
                                                              omega=omega,
                                                              beta=beta,
                                                              mu=mu,
                                                              freeparams=['mu'])
        elif self.MODEL == phydmslib.models.YNGKP_M0:
            e_pw = scipy.asarray([scipy.random.dirichlet([7] * N_NT) for i
                                  in range(3)])
            self.model = phydmslib.models.YNGKP_M0(e_pw, self.nsites)
        else:
            raise ValueError("Invalid MODEL: {0}".format(type(self.MODEL)))

        # make a test tree
        # tree is two sequences separated by a single branch
        # the units are in sub/site
        self.t = 0.04
        newicktree = '(tip1:{0},tip2:{0});'.format(self.t / 2.0)
        self.tree_fname = '_temp.tree'
        with open(self.tree_fname, 'w') as f:
            f.write(newicktree)
        self.tree = Bio.Phylo.read(self.tree_fname, 'newick')

        self.rescaled_tree = copy.deepcopy(self.tree)
        # re-scale the branch lengths
        for node in self.rescaled_tree.get_terminals() + self.rescaled_tree.get_nonterminals():
            if node.branch_length:
                node.branch_length /= self.model.branchScale

        # simulate `Simulator` alignments
        self.simulator_alignments = []
        simulator = phydmslib.simulate.Simulator(self.tree, self.model)
        for i in range(self.nsim):
            self.simulator_alignments.append(simulator.simulate(i))
        # simulate `pyvovle` alignments
        pyvolve_fname = phydmslib.simulate.simulateAlignment(self.model,
                                                             self.tree_fname,
                                                             "{0}_pyvolve".format("test"),
                                                             randomSeed=1,
                                                             nSim=self.nsim)
        self.pyvolve_alignments = [ReadCodonAlignment(fname, checknewickvalid=True)
                                   for fname in pyvolve_fname]
        for fname in pyvolve_fname:
            if os.path.isfile(fname):
                os.remove(fname)

    def test_amino_acid_frequencies(self):
        """Ensure average site-wise amino-acid frequences are similar."""
        simulator_aa_freqs = self.calc_AA_freqs(self.simulator_alignments)
        pyvovlve_aa_freqs = self.calc_AA_freqs(self.pyvolve_alignments)

        # The tolerance is rather lenient because the tree is so small.
        self.assertTrue(scipy.allclose(simulator_aa_freqs,
                                       pyvovlve_aa_freqs, atol=1e-1))

    def test_amino_acid_identity(self):
        """Ensure pair-wise amino-acid identity is similar."""
        simulator_identity = self.calc_AA_identity(self.simulator_alignments)
        pyvolve_identity = self.calc_AA_identity(self.pyvolve_alignments)

        self.assertTrue(scipy.allclose(simulator_identity,
                                       pyvolve_identity, atol=1e-2))

    def test_site_entropy(self):
        """Ensure site-wise entropy is similar."""
        simulator_entropy = self.calc_site_entropy(self.simulator_alignments)
        pyvolve_entropy = self.calc_site_entropy(self.simulator_alignments)

        self.assertTrue(scipy.allclose(simulator_entropy, pyvolve_entropy))

    def calc_site_entropy(self, a):
        """Calculate the average site-wise entropy."""
        # set up
        amino_acids = list(INDEX_TO_AA.values())
        amino_acids.sort()

        site_entropy = [[] for x in range(self.nsites)]
        temp = [self.calc_aa_frequencies(x)[amino_acids].values for x in a]
        for sim in temp:
            for i, site in enumerate(sim):
                site_entropy[i].append(entropy(site))
        site_entropy = [scipy.average(scipy.array(x)) for x in site_entropy]
        return scipy.array(site_entropy)

    def calc_AA_identity(self, a):
        """Calculate average pairwise amino-acid identity."""
        final = []
        for sim in a:
            aa_seq = [self.translate_with_gaps(seq[1]) for seq in sim]
            for seq_pair in itertools.combinations(aa_seq, 2):
                aa_id = [1 if seq_pair[0][i] == seq_pair[1][i] else 0
                         for i in range(len(seq_pair[0]))]
                aa_id = sum(aa_id) / len(aa_id)
                final.append(aa_id)
        final = sum(final) / len(final)
        return final

    def calc_AA_freqs(self, a):
        """Calculate average site-wise amino-acid frequencies."""
        # set up
        amino_acids = list(INDEX_TO_AA.values())
        amino_acids.sort()

        aa_freqs = [[] for x in range(self.nsites)]
        temp = [self.calc_aa_frequencies(x)[amino_acids].values for x in a]
        for sim in temp:
            for i, site in enumerate(sim):
                aa_freqs[i].append(site)
        aa_freqs = scipy.array([scipy.average(scipy.array(x), axis=0)
                                for x in aa_freqs])
        for site in aa_freqs:
            assert scipy.allclose(scipy.sum(site), 1.0)
        return aa_freqs

    def tearDown(self):
        """Remove some files made by `pyvolve`."""
        for f in ['custom_matrix_frequencies.txt']:
            if os.path.isfile(f):
                os.remove(f)
        os.remove(self.tree_fname)

    def translate_with_gaps(self, seq):
        new_seq = []
        for i in range(len(seq) // 3):
            codon = seq[3 * i : 3 * i + 3]
            if codon == '---':
                new_seq.append("-")
            else:
                new_seq.append(CODONSTR_TO_AASTR[codon])
        return("".join(new_seq))

    def calc_aa_frequencies(self, alignment):
        """Calculate amino-acid frequencies from a codon alignment."""
        # Read in the alignnment
        assert scipy.all(scipy.array([len(s[1]) % 3 for s in alignment]) == 0),\
            "At least one sequence in the alignment is not a multiple of 3."
        seqlength = len(alignment[0][1]) // 3
        df = {k: [0 for x in range(seqlength)] for k in list(INDEX_TO_AA.values())}

        # count amino acid frequencies
        for seq in alignment:
            for i in range(seqlength):
                codon = seq[1][3 * i: 3 * i + 3]
                if codon != '---':
                    df[CODONSTR_TO_AASTR[codon]][i] += 1
        df = pd.DataFrame(df)

        # Normalize the dataframe
        assert not scipy.any(df.sum(axis=1) == 0), ("Attempting to normalize a "
                                                    "site by an amino acid count"
                                                    " of zero. Does the alignment"
                                                    " have an all gap column?")
        df = df.div(df.sum(axis=1), axis=0)
        assert scipy.allclose(df.sum(axis=1), 1, atol=0.005)

        # Final formatting
        aa = [x for x in INDEX_TO_AA.values()]
        aa.sort()  # ABC order
        final_cols = ["site"]
        final_cols.extend(aa)
        df["site"] = [x+1 for x in range(len(df))]
        df = df[final_cols]
        return df


class test_compare_YNGKP_M0(test_compare_ExpCM):
    """Tests `pyvolve` and `Simulator` on `YNGKP_M0`."""
    MODEL = phydmslib.models.YNGKP_M0


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
