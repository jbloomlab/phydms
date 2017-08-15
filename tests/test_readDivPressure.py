"""Tests `phydmslib.file_io.readDivPressure`.

Written by Jesse Bloom.
Edited by Sarah Hilton
"""

import os
import unittest
import numpy
import phydmslib.file_io


class test_readPrefs(unittest.TestCase):

    def setUp(self):
        self.divPressurefiles = [os.path.abspath(os.path.join(os.path.dirname
                (__file__), './divpressure_tests/divpressure{0}'.format(suffix)))
                for suffix in ['.csv', '.tsv', '.txt']]
        self.assertTrue(all(map(os.path.isfile, self.divPressurefiles)), "Cannot "
                "find divpressure needed for test.")

    def test_readDivPressure(self):
        """Read divPresures in three different formats."""
        divPressure_as_array = {} # convert to numpy.array so we can use allclose
        for divPressureFile in self.divPressurefiles:
            divPressure = phydmslib.file_io.readDivPressure(divPressureFile)
            sites = sorted(divPressure.keys())
            divPressureList = []
            for r in sites:
            	divPressureList.append(divPressure[r])
            divPressure_as_array[divPressureFile] = numpy.array(divPressureList)
        for (i, f1) in enumerate(self.divPressurefiles):
            divPressure1 = divPressure_as_array[f1]
            for f2 in self.divPressurefiles[i + 1: ]:
                divPressure2 = divPressure_as_array[f2]
                self.assertTrue(divPressure1.shape == divPressure2.shape, "Did not read "
                        "same number of diversifying pressure values from {0} and {1}".format(f1, f2))
                self.assertTrue(numpy.allclose(divPressure1, divPressure2), "Did not get "
                        "same diversifying pressure values from {0} and {1}".format(f1, f2))






if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
