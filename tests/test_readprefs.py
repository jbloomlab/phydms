"""Tests `phydmslib.file_io.readPrefs`.

Written by Jesse Bloom.
"""

import os
import unittest
import numpy
import phydmslib.file_io


class test_readPrefs(unittest.TestCase):

    def setUp(self):
        self.prefsfiles = [os.path.abspath(os.path.join(os.path.dirname
                (__file__), './NP_data/NP_prefs{0}'.format(suffix)))
                for suffix in ['.csv', '.tsv', '-dms_tools_format.txt']]
        self.assertTrue(all(map(os.path.isfile, self.prefsfiles)), "Cannot "
                "find prefsfiles needed for test.")

    def test_readPrefs(self):
        """Read preferences in three different formats."""
        prefs_as_array = {} # convert to numpy.array so we can use allclose
        for prefsfile in self.prefsfiles:
            prefs = phydmslib.file_io.readPrefs(prefsfile)
            sites = sorted(prefs.keys())
            aas = sorted(prefs[sites[0]].keys())
            preflist = []
            for r in sites:
                for aa in aas:
                    preflist.append(prefs[r][aa])
            prefs_as_array[prefsfile] = numpy.array(preflist)
        for (i, f1) in enumerate(self.prefsfiles):
            prefs1 = prefs_as_array[f1]
            for f2 in self.prefsfiles[i + 1: ]:
                prefs2 = prefs_as_array[f2]
                self.assertTrue(prefs1.shape == prefs2.shape, "Did not read "
                        "same number of prefs from {0} and {1}".format(f1, f2))
                self.assertTrue(numpy.allclose(prefs1, prefs2), "Did not get "
                        "same prefs from {0} and {1}".format(f1, f2))

    def test_avgPrefs(self):
        """Tests that the `avgprefs` option works."""
        prefs = phydmslib.file_io.readPrefs(self.prefsfiles[0], avgprefs=True)
        firstsite = True
        for r in list(prefs.keys()):
            if firstsite:
                firstsite = False
                aas = sorted(prefs[r].keys())
                firstsiteprefs = numpy.array([prefs[r][aa] for aa in aas])
            else:
                siteprefs = numpy.array([prefs[r][aa] for aa in aas])
                self.assertTrue(numpy.allclose(firstsiteprefs, siteprefs),
                        "Prefs not same for all sites wth avgprefs")

    def test_minPrefs(self):
        """Tests that the `minprefs` option works."""
        minpref = 0.01
        prefs = phydmslib.file_io.readPrefs(self.prefsfiles[0], minpref=0)
        prefsmin = phydmslib.file_io.readPrefs(self.prefsfiles[0], minpref=minpref)
        sites = sorted(prefs.keys())
        aas = sorted(prefs[sites[0]].keys())
        prefslist = []
        prefsminlist = []
        for r in sites:
            for aa in aas:
                prefslist.append(prefs[r][aa])
                prefsminlist.append(prefsmin[r][aa])
        prefsarray = numpy.array(prefslist)
        prefsminarray = numpy.array(prefsminlist)
        self.assertTrue((prefsarray < minpref).any(), "Not a good test "
                "as all prefs already >= minpref")
        self.assertFalse((prefsminarray < minpref).any(), "Still found "
                "prefs < minpref")



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
