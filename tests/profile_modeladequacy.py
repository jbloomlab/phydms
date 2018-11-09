"""Profiles modeladequacy on a small-ish data set.

Written by Sarah Hilton
"""

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

def profile_pyvolve(alignment_fname, preferences_fname, tree_fname, modelparams_fname, n_sim, n_sites, n_seqs):
    """Main body of script."""
    # inputs
    metrics = ["JensenShannon", "half_sum_abs_diff"]

    # set up
    amino_acids = list(INDEX_TO_AA.values())
    amino_acids.sort()

    # Read in data
    prefs = phydmslib.file_io.readPrefs(preferences_fname)
    prefs = [prefs[r] for r in sorted(prefs.keys())]
    alignment = phydmslib.file_io.ReadCodonAlignment(alignment_fname,
                                                     checknewickvalid=True)
    alignment_freqs = phydmslib.modeladequacy.calc_aa_frequencies(alignment)
    alignment_freqs["alignment"] = (alignment_freqs[amino_acids]
                                    .apply(lambda r: tuple(r), axis=1)
                                    .apply(scipy.array))
    alignment_freqs = scipy.array(alignment_freqs["alignment"].tolist())

    # build model
    model = phydmslib.modeladequacy.make_expcm(modelparams_fname, prefs)
    # stationary state frequencies
    ss_freqs = phydmslib.modeladequacy.calc_stationary_state_freqs(model)

    # simulate sequences
    sim_fname = phydmslib.simulate.simulateAlignment(model,  # expcm
                                                     tree_fname,  # tree
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

    # process simulations
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
    pd.concat(final).to_csv("_modeladequacy_pvalues_nsites{0}_nseqs{1}_pyvovle.csv".format(n_sites, n_seqs), index=False)

    for fname in glob.glob("_modeladequacy_results_*"):
        os.remove(fname)

def profile_simulator(alignment_fname, preferences_fname, tree_fname, modelparams_fname, n_sim, n_sites, n_seqs):
    """Main body of script."""
    # inputs
    metrics = ["JensenShannon", "half_sum_abs_diff"]

    # set up
    amino_acids = list(INDEX_TO_AA.values())
    amino_acids.sort()

    # Read in data
    prefs = phydmslib.file_io.readPrefs(preferences_fname)
    prefs = [prefs[r] for r in sorted(prefs.keys())]
    alignment = phydmslib.file_io.ReadCodonAlignment(alignment_fname,
                                                     checknewickvalid=True)
    alignment_freqs = phydmslib.modeladequacy.calc_aa_frequencies(alignment)
    alignment_freqs["alignment"] = (alignment_freqs[amino_acids]
                                    .apply(lambda r: tuple(r), axis=1)
                                    .apply(scipy.array))
    alignment_freqs = scipy.array(alignment_freqs["alignment"].tolist())

    # build model
    model = phydmslib.modeladequacy.make_expcm(modelparams_fname, prefs)
    # stationary state frequencies
    ss_freqs = phydmslib.modeladequacy.calc_stationary_state_freqs(model)

    # read in tree
    tree = Bio.Phylo.read(tree_fname, 'newick')
    for node in tree.get_terminals() + tree.get_nonterminals():
        if (node.branch_length is None) and (node == tree.root):
            node.branch_length = 1e-06

    # simulate the sequences
    simulated_sequences = []
    simulator = phydmslib.simulate.Simulator(tree, model)
    for fname in range(n_sim):
        simulated_sequences.append(simulator.simulate())
    # simulated sequences is a list of lists
    # There is one sequence list for each simulation replicate
    # each sequence list is a list of tuples, one tuple per sequence

    # need to process simulations on a site by site basis
    simulations = [[] for x in range(model.nsites)]
    for sim in simulated_sequences: # for each simulation replicate
        sim_freqs = phydmslib.modeladequacy.calc_aa_frequencies(sim) # each row is a site and each column is an amino acid
        sim_freqs = sim_freqs[amino_acids].values # each row is a site and each column is an amino acid
        for site in range(model.nsites): # for each site, put the amino acid frequency array into a new list
            simulations[site].append(sim_freqs[site])

    # process simulations
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
    pd.concat(final).to_csv("_modeladequacy_pvalues_nsites{0}_nseqs{1}_simulator.csv".format(n_sites, n_seqs), index=False)

if __name__ == '__main__':
    for nsite in [10]:
        for nseq in [34]:
            alignment = "modeladequacy_tests/HA_short_nsites{0}_nseqs{1}.fasta".format(nsite, nseq)
            prefs = "modeladequacy_tests/HA_short_prefs_nsites{0}.csv".format(nsite)
            tree = "modeladequacy_tests/HA_short_tree_nsites{0}_nseqs{1}.newick".format(nsite, nseq)
            mp = "modeladequacy_tests/phydms_nsites{0}_nseqs{1}/modelparams.txt".format(nsite, nseq)
            # first with `pyvolve`
            statsfile = 'cprofile_modeladequacy_nsites{0}_nseqs{1}_pyvolve'.format(nsite, nseq)
            pstatsfile = "pstats_modeladequacy_nsites{0}_nseqs{1}_pyvolve".format(nsite, nseq)
            cProfile.run('profile_pyvolve(alignment, prefs, tree, mp, 1000, nsite, nseq)', statsfile)
            with open(pstatsfile, 'w') as stream:
                stats = pstats.Stats(statsfile, stream=stream)
                for t in ['cumtime', 'tottime']:
                    print(f'\n{t}:')
                    # stats.strip_dirs().sort_stats(t).print_stats(100)
                    stats.sort_stats(t).print_stats(100)

            # second with `simulator`
            statsfile = 'cprofile_modeladequacy_nsites{0}_nseqs{1}_simulator'.format(nsite, nseq)
            pstatsfile = "pstats_modeladequacy_nsites{0}_nseqs{1}_simulator".format(nsite, nseq)
            cProfile.run('profile_simulator(alignment, prefs, tree, mp, 1000, nsite, nseq)', statsfile)
            with open(pstatsfile, 'w') as stream:
                stats = pstats.Stats(statsfile, stream=stream)
                for t in ['cumtime', 'tottime']:
                    print(f'\n{t}:')
                    # stats.strip_dirs().sort_stats(t).print_stats(100)
                    stats.sort_stats(t).print_stats(100)
