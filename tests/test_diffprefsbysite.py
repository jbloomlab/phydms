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
import scipy.stats
import phydmslib.file_io
import phydmslib.models
import phydmslib.simulate
from phydmslib.constants import *
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pyvolve


class test_OmegaBySiteExpCM(unittest.TestCase):
    """Tests ``--omegabysite`` to ``phydms`` for `ExpCM`."""
    gammaomega_arg = []

    def setUp(self):
        random.seed(1)
        scipy.random.seed(1)
        self.tree = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_tree.newick'))
        self.alignment = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_alignment.fasta'))
        self.prefs = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_prefs.tsv'))
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
        self.outdir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './diffprefsbysite_test_results/'))
        if self.gammaomega_arg:
            self.modelname = 'ExpCM_gammaomega'
        else:
            self.modelname = 'ExpCM'
        self.modelarg = 'ExpCM_{0}'.format(self.prefs)

        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)

    def test_OnSimulatedData(self):
        """Run ``phydms`` on the simulated data."""
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

        prefsbymethod = {}
        for fitprefsmethod in ['1', '2']:
            outprefix = simulateprefix + '_fitprefsmethod{0}'.format(
                    fitprefsmethod)
            subprocess.check_call(['phydms', simulatedalignment, self.tree,
                    self.modelarg, outprefix, '--diffprefsbysite',
                    '--brlen', 'scale', '--ncpus', '-1', '--diffprefsprior',
                    'invquadratic,150,0.5'] + self.gammaomega_arg +
                    ['--fitprefsmethod', fitprefsmethod])
            diffprefsbysitefile = outprefix + '_diffprefsbysite.txt'
            aas = ['dpi_{0}'.format(INDEX_TO_AA[a]) for a in range(N_AA)]
            diffprefs = pandas.read_csv(diffprefsbysitefile, sep='\t',
                    comment='#')
            diffprefs['total'] = diffprefs[aas].abs().sum(axis=1)
            for (site, a) in self.targetaas.items():
                siteentry = diffprefs[diffprefs['site'] == site]
                self.assertTrue(len(siteentry) == 1, str(len(siteentry)))
                self.assertTrue((siteentry['dpi_{0}'.format(a)] > 0).all())

            prefsbymethod[fitprefsmethod] = diffprefs

        for (i, (method1, prefs1)) in enumerate(sorted(prefsbymethod.items())):
            total1 = prefs1['total'].values
            for (method2, prefs2) in sorted(prefsbymethod.items())[i + 1 : ]:
                total2 = prefs2['total'].values
                (r, p) = scipy.stats.pearsonr(total1, total2)
                plt.scatter(total1, total2)
                plt.xlabel('fitprefsmethod{0}'.format(method1))
                plt.ylabel('fitprefsmethod{0}'.format(method2))
                plotfile = os.path.join(self.outdir, '{0}_vs_{1}.pdf'.format(
                        method1, method2))
                plt.savefig(plotfile)
                self.assertTrue(r > 0.98, "Low correlation between "
                        "fitprefsmethods: {0}\nSee {1}".format(r, plotfile))

        for f in ["custom_matrix_frequencies.txt"]:
            if os.path.isfile(f):
                os.remove(f)



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
