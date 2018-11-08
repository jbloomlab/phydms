"""Functions for simulating alignments.

Simulations can be performed by using `pyvolve` partitions or `phydms` models.

Written by Jesse Bloom and Sarah Hilton.
"""

import os
import sys
import scipy
import math
import phydmslib.models
from phydmslib.constants import *
import pyvolve
from tempfile import mkstemp
import random
import Bio.Phylo
import numpy as np
import copy


class Simulator(object):
    """Uses model and tree to simulate an alignment.

    This class uses the `phydms` models to simulate an alignment.
    Most commonly, you will initiate a simulator object, simulate
    an alignment, and then output the final alignment.
    > simulator = phydmslib.simulate.simulator(tree, model)
    > simulator.simulate()
    > simulator.output_alignment(fname)

    See `__init__` for how to initialize a `simulator`.

    Attributes:
        `tree` (instance of `Bio.Phylo.BaseTree.Tree` derived class)
            Phylogenetic tree. Branch lengths are in units of codon
            substitutions per site for the current `model`.
        `model` (instance of `phydmslib.models.Model` derived class)
            Specifies the substitution model for `nsites` codon sites.
            Only accepts `ExpCM_empirical_phi`, `ExpCM`, `YNGKP_M0`
        `nsites` (int)
            Number of codon sites.
        `root` (instance of `Bio.Phylo.BaseTree.Tree.clade` derived class)
            The root of `tree`
        `internalnode` (list)
            Names of the internal nodes in `tree`.
        `terminalnode` (list)
            Names of the terminal nodes in `tree`.
        `alignment` (dictionary)
            Sequence arrays for each node in the `tree`. The dictionary is
            keyed by the node name and the value is a list of `nsites`
            `scipy.arrays` of length `N_CODON`. The codon identity at the site
            is a 1 and all other elements of the array are zeroes.
        `simulated_alignment` (list)
            Final codon alignment. Tuples of (seq_id, seq) for each node in
            the `tree`.
    """
    def __init__(self, tree, model, branchScale=None):
        """Initialize a `simulator` object.

        Args:
            `tree`, `model`, `model`
                Attributes of same name described in class doc string.
                Note that we make copies of `tree`, and `model`, so
                the calling objects are not modified during simulation.
            `branchScale` (`None` or float > 0)
                The branch lengths in the input tree are assumed
                to be in units of substitutions per site with
                the scaling defined by `model.branchScale`. If
                the scaling should instead be defined by some
                other value of `branchScale`, indicate by
                setting to that value rather than `None`. This
                is useful if tree was inferred on models across
                many sites but you are now just analyzing an
                individual one. Note that this option applies
                only to the input tree, not the output one.
        """
        if isinstance(model, phydmslib.models.ExpCM_empirical_phi) or isinstance(model, phydmslib.models.YNGKP_M0) or isinstance(model, phydmslib.models.ExpCM):
            self.model = copy.deepcopy(model)
            self.nsites = self.model.nsites
        else:
            raise ValueError("model is not a valid model")

        # set up dictionary to hold sequence arrays
        self.alignment = {}
        self.simulated_alignment = []

        # tree - not changing the units. Assuming they are correct
        assert isinstance(tree, Bio.Phylo.BaseTree.Tree), "invalid tree"
        self._tree = copy.deepcopy(tree)

        # For internal storage, branch lengths are in model units rather
        # than codon substitutions per site. So here adjust them from
        # codon substitutions per site to model units.
        self.root = False
        self.internalnode = []
        self.terminalnode = []
        assert isinstance(tree, Bio.Phylo.BaseTree.Tree), "invalid tree"
        self._tree = copy.deepcopy(tree)
        if branchScale == None:
            branchScale = self.model.branchScale
        else:
            assert isinstance(branchScale, float) and branchScale > 0
        for node in self._tree.find_clades():
            if node != self._tree.root:
                node.branch_length /= branchScale
                self.alignment[node] = []
                if len(node.clades) == 0:
                    self.terminalnode.append(node.name)
                else:
                    self.internalnode.append(node.name)
            else:
                assert not self.root, "More than one root!"
                self.root = node
                self.alignment[node] = []

    def simulate(self):
        """Simulate an alignment.

        Calls functions to simulate an alignment and the
        processes the simulated alignment.

        Sequences for a given site are stored as a `scipy.array` of
        length 61. The codon sequence for that site is recorded as a
        1 at the codon index in this array. All other elements are 0s.


        Two private functions:
            `_evolve_branch`: generate a new sequence. Called by `_traverse_tree`
            `_traverse_tree`: recursively walks from root to tip of tree.
        """

        def _evolve_branch(parent, child, alignment):
            """Generates new sequence given a parent and a child node."""
            branch_length = parent.distance(child)
            exp_M = self.model.M(branch_length)
            new_seq = []
            for r in range(self.model.nsites):
                site_distribution = alignment[parent][r].dot(exp_M[r])  # parent's sequence array times the transition matrix
                codon = scipy.random.choice(N_CODON, 1, p=site_distribution)[0]  # choose from the new codon distribution for the site
                # create new sequence array. 1 indicates codon sequence
                site_seq = scipy.zeros(N_CODON)
                site_seq[codon] = 1
                assert len(scipy.where(site_seq==1)[0]) == 1
                new_seq.append(site_seq)
            alignment[child] = scipy.array(new_seq)  # array of size `nsites`, each element array of 61
            return alignment

        def _traverse_tree(parent, child, alignment):
            """Walks along a branch of the tree.
            Calls the `_evolve_branch` function for internal nodes.
            """
            if parent == None:  # at the root, need to generate a sequence
                root_seq = []
                for r in range(self.model.nsites):
                    ss = self.model.stationarystate[r]
                    codon = scipy.random.choice(N_CODON, 1, p=ss)[0]  # pull a codon from the stationary state
                    # create sequence array. 1 indicates codon sequence
                    site_seq = scipy.zeros(N_CODON)
                    site_seq[codon] = 1
                    assert len(scipy.where(site_seq==1)[0]) == 1
                    root_seq.append(site_seq)
                alignment[child] = scipy.array(root_seq)  # array of size `nsites`, each element array of 61
            else:  # internal branch, need to evolve along the branch
                alignment = _evolve_branch(parent, child, alignment)

            # Did we just simulate to a terminal node?
            # If not, call the `_traverse_tree` function again
            if not child.is_terminal():
                for gchild in child.clades:
                    alignment = _traverse_tree(child, gchild, alignment)
            return alignment

        # beginning of `simulate` function
        alignment = {}
        alignment = _traverse_tree(None, self.root, alignment)  # simulate the sequences

        # reformat the simulated alignment
        # turn the sequenc arrays into codon sequnces
        # format alignment as a list of tuples
        simulated_alignment = []
        for key in alignment:
            seq = alignment[key]
            final = []
            for site in seq:
                codon = scipy.where(site==1)[0]
                assert len(codon) == 1, "There is more than one codon at a site!"
                nt = INDEX_TO_CODON[codon[0]]
                final.append(nt)
            final = "".join(final)
            simulated_alignment.append((key.name, final))
        return simulated_alignment

    def output_alignment(self, simulated_alignment, fname):
        """Outputs simulated alignment as fasta file.

        Args:
            `fname`
                Name of outputput fasta file.
        """
        with open(fname, "w") as f:
            for seq in simulated_alignment:
                if seq[0] in self.terminalnode:
                    f.write(">{0}\n{1}\n".format(seq[0],seq[1]))

# simulatoins with `pyvolve`
def pyvolvePartitions(model, divselection=None, rateMatrixPrefix=""):
    """Get list of `pyvolve` partitions for `model`.

    Args:
        `model` (`phydmslib.models.Models` object)
            The model used for the simulations. Currently only
            certain `Models` are supported (e.g., `YNGKP`,
            `ExpCM`)
        `divselection` (`None` or 2-tuple `(divomega, divsites)`)
            Set this option if you want to simulate a subset of sites
            as under diversifying selection (e.g., an `omega` different
            than that used by `model`. In this case, `divomega` is
            the omega for this subset of sites, and `divsites` is a list
            of the sites in 1, 2, ... numbering.

    Returns:
        `partitions` (`list` of `pyvolve.Partition` objects)
            Can be fed into `pyvolve.Evolver` to simulate evolution.

    """
    codons = pyvolve.genetics.Genetics().codons
    codon_dict = pyvolve.genetics.Genetics().codon_dict
    pyrims = pyvolve.genetics.Genetics().pyrims
    purines = pyvolve.genetics.Genetics().purines

    if divselection:
        (divomega, divsites) = divselection
    else:
        divsites = []

    assert all([1 <= r <= model.nsites for r in divsites])

    partitions = []
    for r in range(model.nsites):
        matrix = scipy.zeros((len(codons), len(codons)), dtype='float')
        for (xi, x) in enumerate(codons):
            for (yi, y) in enumerate(codons):
                ntdiffs = [(x[j], y[j]) for j in range(3) if x[j] != y[j]]
                if len(ntdiffs) == 1:
                    (xnt, ynt) = ntdiffs[0]
                    qxy = 1.0
                    if (xnt in purines) == (ynt in purines):
                        qxy *= model.kappa
                    (xaa, yaa) = (codon_dict[x], codon_dict[y])
                    fxy = 1.0
                    if xaa != yaa:
                        if type(model) == phydmslib.models.ExpCM_empirical_phi_divpressure:
                            fxy *= model.omega * (1 + model.omega2 * model.deltar[r])
                        elif r + 1 in divsites:
                            fxy *= divomega
                        else:
                            fxy *= model.omega
                    if type(model) in [phydmslib.models.ExpCM,
                                       phydmslib.models.ExpCM_empirical_phi,
                                       phydmslib.models.ExpCM_empirical_phi_divpressure]:
                        qxy *= model.phi[NT_TO_INDEX[ynt]]
                        pix = model.pi[r][AA_TO_INDEX[xaa]]**model.beta
                        piy = model.pi[r][AA_TO_INDEX[yaa]]**model.beta
                        if abs(pix - piy) > ALMOST_ZERO:
                            fxy *= math.log(piy / pix) / (1.0 - pix / piy)
                    elif type(model) == phydmslib.models.YNGKP_M0:
                        for p in range(3):
                            qxy *= model.phi[p][NT_TO_INDEX[y[p]]]
                    else:
                        raise ValueError("Can't handle model type {0}".format(
                                type(model)))
                    matrix[xi][yi] = model.mu * qxy * fxy
            matrix[xi][xi] = -matrix[xi].sum()

        # create model in way that captures `pyvovle` print statements
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        try:
            custom_matrix_fname = "{0}custom_matrix_frequencies.txt"\
                                  .format(rateMatrixPrefix)
            m = pyvolve.Model("custom",
                              {"matrix": matrix},
                              save_custom_frequencies=custom_matrix_fname)
        finally:
            sys.stdout.close()
            sys.stdout = old_stdout
        partitions.append(pyvolve.Partition(models=m, size=1))

    return partitions


def simulateAlignment(model, treeFile, alignmentPrefix,
                      randomSeed=False, nSim=1):
    """
    Simulate an alignment given a model and tree (units = subs/site).

    Simulations done using `pyvolve`.

    Args:
        `model` (`phydmslib.models.Models` object)
            The model used for the simulations. Only
            models that can be passed to `pyvolve.Partitions`
            are supported.
        `treeFile` (str)
            Name of newick file used to simulate the sequences.
            The branch lengths should be in substitutions per site,
            which is the default units for all `phydms` outputs.
        `alignmentPrefix` (str)
            Prefix for the files created by `pyvolve`.
        `randomSeed` (int)
            The integer used to seed the random number generator.
        `nSim` (int)
            The number of replicate simlulations to perform. If no
            number is specified than only one simulated alignment will
            be created.

    Returns:
        `createdAlignments` (`str` or `list` of simulation file names)
            A list of the output files created if `nSim` > 1 or the
            single output file if `nSim` = 1.

    The result of this function is a simulated FASTA alignment. If `nSim`
    is set to 1 or not specified then the simulated alignment will be named
     `'{alignmentPrefix}_simulatedalignment.fasta'`. Otherwise, there will be
     an alignment named `'{alignmentPrefix}_{rep}_simulatedalignment.fasta'`
     for `rep` in `range(nSim)`.
    """
    if randomSeed is False:
        pass
    else:
        random.seed(randomSeed)
        np.random.seed(randomSeed)

    # Transform the branch lengths by dividing by the model `branchScale`
    tree = Bio.Phylo.read(treeFile, 'newick')
    for node in tree.get_terminals() + tree.get_nonterminals():
        if (node.branch_length is None) and (node == tree.root):
            node.branch_length = 1e-06
        else:
            node.branch_length /= model.branchScale
    fd, temp_path = mkstemp()
    Bio.Phylo.write(tree, temp_path, 'newick')
    os.close(fd)
    pyvolve_tree = pyvolve.read_tree(file=temp_path)
    os.remove(temp_path)

    # Make the `pyvolve` partition
    partitions = pyvolvePartitions(model, rateMatrixPrefix=alignmentPrefix)

    createdAlignments = []
    for rep in range(nSim):
        # Simulate the alignment
        alignment = '{0}_{1}_simulatedalignment.fasta'\
                    .format(alignmentPrefix, rep)
        info = '_temp_{0}info.txt'.format(alignmentPrefix)
        rates = '_temp_{0}_ratefile.txt'.format(alignmentPrefix)
        custom_matrix = '{0}custom_matrix_frequencies.txt'\
                        .format(alignmentPrefix)
        evolver = pyvolve.Evolver(partitions=partitions, tree=pyvolve_tree)
        evolver(seqfile=alignment, infofile=info, ratefile=rates)
        # Clean up extraneous `pyvovle` simulations
        for f in [rates, info, custom_matrix]:
            if os.path.isfile(f):
                os.remove(f)
        assert os.path.isfile(alignment)
        createdAlignments.append(alignment)

    # Make sure all of the expected output files exist
    if nSim == 1:
        # Re-name files if only one replicate is specified
        createdAlignments = ['{0}_simulatedalignment.fasta'\
                             .format(alignmentPrefix)]
        os.rename(alignment, createdAlignments[0])
        assert os.path.isfile(createdAlignments[0]),\
            "Failed to created file {0}".format(createdAlignments[0])
    else:
        for f in createdAlignments:
            assert os.path.isfile(f), "Failed to created file {0}".format(f)
    return createdAlignments


if __name__ == '__main__':
    import doctest
    doctest.testmod()
