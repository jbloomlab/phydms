"""Tests ``--omegabysite`` option to  ``phydms`` on simulate data.

Written by Jesse Bloom.
"""

import os
import shutil
import unittest
import subprocess
import random
import pandas
import phydmslib.file_io
import phydmslib.models
import phydmslib.simulate
from phydmslib.constants import *
import pyvolve


class test_OmegaBySiteExpCM(unittest.TestCase):
    """Tests ``--omegabysite`` to ``phydms`` for `ExpCM`."""

    def setUp(self):
        self.tree = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_tree.newick'))
        self.alignment = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_alignment.fasta'))
        self.prefs = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_prefs.tsv'))
        self.nsites = len(phydmslib.file_io.ReadCodonAlignment(self.alignment,
                True)[0][1]) // 3
        self.initializeModel()
        self.outdir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './omegabysite_test_results/'))
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)

    def initializeModel(self):
        prefs = phydmslib.file_io.readPrefs(self.prefs, minpref=0.005)
        prefs = [prefs[r] for r in sorted(list(prefs.keys()))]
        # Using beta < 1 partially flattens prefs in simulation
        # Use mu < 1 to get branch lengths about right
        self.model = phydmslib.models.ExpCM(prefs, beta=0.5, mu=0.5)
        self.modelname = 'ExpCM'
        self.modelarg = 'ExpCM_{0}'.format(self.prefs)

    def test_OnSimulatedData(self):
        random.seed(1)
        divpressuresites = random.sample(range(self.nsites), 5)
        partitions = phydmslib.simulate.pyvolvePartitions(self.model,
                (200.0, divpressuresites))
        evolver = pyvolve.Evolver(partitions=partitions,
                tree=pyvolve.read_tree(file=self.tree))
        simulateprefix = os.path.join(self.outdir, self.modelname)
        simulatedalignment = simulateprefix + '_simulatedalignment.fasta'
        info = simulateprefix + '_temp_info.txt'
        rates = simulateprefix + '_temp_ratefile.txt'
        evolver(seqfile=simulatedalignment, infofile=info, ratefile=rates)
        subprocess.check_call(['phydms', simulatedalignment, self.tree,
                self.modelarg, simulateprefix, '--omegabysite',
                '--brlen', 'scale'])
        omegabysitefile = simulateprefix + '_omegabysite.txt'
        omegas = pandas.read_csv(omegabysitefile, sep='\t', comment='#')
        divpressureomegas = omegas[omegas['site'].isin(divpressuresites)]
        self.assertTrue(len(divpressureomegas) == len(divpressuresites))
        self.assertTrue((divpressureomegas['omega'].values > 2).all(),
                "Not all divpressure sites have omega > 2:\n{0}".format(
                divpressureomegas))
        self.assertTrue((divpressureomegas['P'].values < 0.08).all(),
                "Not all divpressure sites have P < 0.08:\n{0}".format(
                divpressureomegas))
        nspurious = len(omegas[(omegas['omega'] > 2) & (omegas['P'] < 0.05)
                & (~omegas['site'].isin(divpressuresites))])
        self.assertTrue(nspurious <= 1, "{0} spurious sites".format(nspurious))

        for f in ["custom_matrix_frequencies.txt"]:
            if os.path.isfile(f):
                os.remove(f)


class test_OmegaBySiteYNGKP(test_OmegaBySiteExpCM):
    """Tests ``--omegabysite`` to ``phydms`` for `YNGKP_M0`."""

    def initializeModel(self):
        e_pw = scipy.full((3, N_NT), 1.0 / N_NT, dtype='float')
        # mu > 1 leads to longer branches in simulation
        self.model = phydmslib.models.YNGKP_M0(e_pw, self.nsites, mu=4.0)
        self.modelname = 'YNGKP'
        self.modelarg = 'YNGKP_M0'


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
