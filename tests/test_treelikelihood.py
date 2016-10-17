"""Tests `phydmslib.treelikelihood.TreeLikelihood`.

Written by Jesse Bloom.
"""

import os
import sys
import re
import unittest
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import Bio.Phylo


class test_TreeLikelihood(unittest.TestCase):
    """Tests `phydmslib.treelikelihood.TreeLikelihood` class."""

    def setUp(self):
        """Set up parameters for test."""

        # define tree and write image to a file
        self.newick = ('(($node_1=CAA$:0.1,$node_2=CAG$:0.15)$node_4=x$:0.15,'
                       '$node_3=GAA$:0.25)$node_5=y$:0.02;')
        tempfile = '_temp.tree'
        with open(tempfile, 'w') as f:
            f.write(self.newick)
        self.tree = Bio.Phylo.read(tempfile, 'newick')
        os.remove(tempfile)
        branch_labels = {} # branch annotations
        cladebyname = dict([(clade.name, clade) for clade in 
                self.tree.find_clades()])
        for (name, brlen) in re.findall(
                '(?P<name>\$node_\d\=[A-z]+\$)\:(?P<brlen>\d+\.\d+)', 
                self.newick):
            if name != self.tree.root.name:
                i = name.split('=')[0][-1] # node number
                branch_labels[cladebyname[name]] = "$t_{0}={1}$".format(
                        i, brlen)
        matplotlib.rc('text', usetex=True)
        Bio.Phylo.draw(self.tree, do_show=False, branch_labels=branch_labels)
        plt.axis('off')
        plt.savefig('test_treelikelihood_image.pdf')

    def test_TreeLikelihood(self):
        pass


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
