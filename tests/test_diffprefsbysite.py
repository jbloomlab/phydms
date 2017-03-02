"""Tests ``--diffprefsbysite`` option to  ``phydms`` on simulated data.

Written by Jesse Bloom.
"""

import os
import shutil
import unittest
import subprocess
import random
import pandas
import scipy
import phydmslib.file_io
import phydmslib.models
import phydmslib.simulate
from phydmslib.constants import *
import pyvolve


class test_OmegaBySiteExpCM(unittest.TestCase):
    """Tests ``--omegabysite`` to ``phydms`` for `ExpCM`."""

    def setUp(self):
        random.seed(1)
        scipy.random.seed(1)
        self.tree = './NP_data/NP_tree.newick'
        self.alignment = './NP_data/NP_alignment.fasta'
        self.prefs = './NP_data/NP_prefs.tsv'
        self.nsites = len(phydmslib.file_io.ReadCodonAlignment(self.alignment, 
                True)[0][1]) // 3
        prefs = phydmslib.file_io.readPrefs(self.prefs, minpref=0.005)
        aas = [INDEX_TO_AA[a] for a in range(N_AA)]
        self.shuffledsites = random.sample(sorted(prefs.keys()), 10)
        self.targetaas = dict([(r, random.choice(aas)) for r in self.shuffledsites])
        prefs = [prefs[r] for r in sorted(list(prefs.keys()))]
        hipref = 0.9
        lowpref = (1.0 - hipref) / (N_AA - 1)
        for r in self.shuffledsites:
            rprefs = [lowpref] * N_AA
            rprefs[AA_TO_INDEX[self.targetaas[r]]] = hipref
            prefs[r - 1] = dict(zip(aas, rprefs))
        self.model = phydmslib.models.ExpCM(prefs, beta=1.5)
        self.outdir = './diffprefsbysite_test_results/'
        self.modelname = 'ExpCM'
        self.modelarg = 'ExpCM_{0}'.format(self.prefs)

        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)

    def test_OnSimulatedData(self):
        random.seed(1)
        scipy.random.seed(1)
        partitions = phydmslib.simulate.pyvolvePartitions(self.model)
        evolver = pyvolve.Evolver(partitions=partitions, 
                tree=pyvolve.read_tree(file=self.tree))
        simulateprefix = os.path.join(self.outdir, self.modelname)
        simulatedalignment = simulateprefix + '_simulatedalignment.fasta'
        info = simulateprefix + '_temp_info.txt'
        rates = simulateprefix + '_temp_ratefile.txt'
        evolver(seqfile=simulatedalignment, infofile=info, ratefile=rates)
        subprocess.check_call(['phydms', simulatedalignment, self.tree,
                self.modelarg, simulateprefix, '--diffprefsbysite', 
                '--brlen', 'scale', '--ncpus', '-1', '--diffprefsprior',
                'invquadratic,50,0.25'])
        diffprefsbysitefile = simulateprefix + '_diffprefsbysite.txt'
        aas = [INDEX_TO_AA[a] for a in range(N_AA)]
        diffprefs = pandas.read_csv(diffprefsbysitefile, sep='\t', comment='#',
                header=None, names=(['site'] + aas))
        diffprefs['total'] = diffprefs[aas].sum(axis=1)
        diffprefs.sort_values('total', ascending=False, inplace=True)
        print(diffprefs)
        print(self.shuffledsites)
        print(self.targetaas)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
