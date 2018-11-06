"""Runs the model adequacy protocol.

This test examines the functionality of the model adequcy protocol and checks
the results against the output in `expected_modeladequacy_results`.

Written by Sarah Hilton.
"""

import unittest
import phydmslib.simulate
import phydmslib.file_io
import phydmslib.modeladequacy
from phydmslib.constants import *
import os
from statsmodels.sandbox.stats.multicomp import multipletests
import scipy
import numpy as np
import itertools
import pandas as pd
import cProfile
import pstats
import glob


class test_modeladequacy_largeTree(unittest.TestCase):
    """Tests the model adequacy protocol on a large tree."""
    # run parameters
    ALIGNMENT = "modeladequacy_tests/HA_short_nsites10_nseqs34.fasta"
    PREFERNCES = "modeladequacy_tests/HA_short_prefs_nsites10.csv"
    TREE = "modeladequacy_tests/HA_short_tree_nsites10_nseqs34.newick"
    MODELPARAMS = "modeladequacy_tests/phydms_nsites10_nseqs34/modelparams.txt"
    EXPECTED = "expected_modeladequacy_results/_modeladequacy_pvalues_nsites10_nseqs34.csv"

    def test_modeladequacy(self):
        """Tests model adequacy on small HA data."""
        # inputs
        n_sim = 10
        metrics = ["JensenShannon", "half_sum_abs_diff"]

        # set up
        amino_acids = list(INDEX_TO_AA.values())
        amino_acids.sort()

        # Read in data
        prefs = phydmslib.file_io.readPrefs(self.PREFERNCES)
        prefs = [prefs[r] for r in sorted(prefs.keys())]
        alignment = phydmslib.file_io.ReadCodonAlignment(self.ALIGNMENT,
                                                         checknewickvalid=True)
        alignment_freqs = phydmslib.modeladequacy.calc_aa_frequencies(alignment)
        alignment_freqs["alignment"] = (alignment_freqs[amino_acids]
                                        .apply(lambda r: tuple(r), axis=1)
                                        .apply(scipy.array))
        alignment_freqs = scipy.array(alignment_freqs["alignment"].tolist())

        # build model
        model = phydmslib.modeladequacy.make_expcm(self.MODELPARAMS, prefs)
        # stationary state frequencies
        ss_freqs = phydmslib.modeladequacy.calc_stationary_state_freqs(model)

        # simulate sequences
        sim_fname = phydmslib.simulate.simulateAlignment(model,  # expcm
                                                         self.TREE,  # tree
                                                         "_modeladequacy_results",  # alignment name
                                                         0,  # seed
                                                         n_sim)  # number of rep
        simulations = [[] for x in range(model.nsites)]
        for sim in sim_fname:
            sim_freqs = phydmslib.file_io.ReadCodonAlignment(sim, checknewickvalid=True)
            sim_freqs = phydmslib.modeladequacy.calc_aa_frequencies(sim_freqs)
            sim_freqs = sim_freqs[amino_acids].values
            for site in range(model.nsites):
                simulations[site].append(sim_freqs[site])

        # process simulations, calculate distances, and pvalues
        np.random.seed(0)
        final = []
        for site in range(model.nsites):
            sims = simulations[site]
            ss = ss_freqs[site]
            natural = alignment_freqs[site]
            for metric in metrics:
                # calc simulation distances
                sim = scipy.array([phydmslib.modeladequacy.prefDistance(sim, ss, metric)
                                   for sim in sims])
                # # calc distances natural
                natural_distance = phydmslib.modeladequacy.prefDistance(natural, ss,
                                                                        metric)
                # calc pvalues
                greater = scipy.sum(scipy.greater(sim, natural_distance))
                tie_breaker = scipy.sum(scipy.equal(sim, natural_distance))
                if tie_breaker >= 1:
                    tie_breaker = np.random.randint(tie_breaker, size=1)[0]
                pvalue = (greater + tie_breaker + 1) / (len(sim) + 1)
                assert 0 <= pvalue <= 1.0, "pvalue is > 1.0 or < 0.0"
                final.append((site, pvalue, metric))

        # format pvalues and calculate qvalues
        df = pd.DataFrame(final, columns=['site', 'pvalue', 'metric'])
        final = []
        for name, group in df.groupby("metric"):
            group = group.sort_values(by="pvalue", ascending=True)
            group["qvalue"] = multipletests(group["pvalue"].tolist(),
                                            method='fdr_bh')[1]
            final.append(group)
        final = pd.concat(final).sort_values(by=["site", "metric"]).drop(['metric'], axis=1)

        # compare results
        expected = pd.read_csv(self.EXPECTED).sort_values(by=["site", "metric"]).drop(['metric'], axis=1)
        self.assertTrue(scipy.allclose(final, expected))

        # remove files
        for fname in glob.glob("_modeladequacy_results_*"):
            os.remove(fname)

class test_modeladequacy_smallTree(test_modeladequacy_largeTree):
    """Tests model adequacy on a small tree."""
    # run parameters
    ALIGNMENT = "modeladequacy_tests/HA_short_nsites10_nseqs5.fasta"
    PREFERNCES = "modeladequacy_tests/HA_short_prefs_nsites10.csv"
    TREE = "modeladequacy_tests/HA_short_tree_nsites10_nseqs5.newick"
    MODELPARAMS = "modeladequacy_tests/phydms_nsites10_nseqs5/modelparams.txt"
    EXPECTED = "expected_modeladequacy_results/_modeladequacy_pvalues_nsites10_nseqs5.csv"


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
