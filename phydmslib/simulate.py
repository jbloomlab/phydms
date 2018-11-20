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
    """Uses a model and a tree to simulate an alignment.

    This class uses the `phydms` models to simulate an alignment.
    Most commonly, you will initiate a `Simulator` object and then
    simulate an aligment using `Simulator.simulate`.

    See `__init__` for how to initialize a `Simulator`.

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

    """

    def __init__(self, tree, model, branchScale=None):
        """Initialize a `simulator` object.

        Args:
            `tree`, `model`
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
                individual one.

        """
        if (isinstance(model, phydmslib.models.ExpCM_empirical_phi)
               or isinstance(model, phydmslib.models.YNGKP_M0)
               or isinstance(model, phydmslib.models.ExpCM)):
            self._model = copy.deepcopy(model)
            self.nsites = self._model.nsites
        else:
            raise ValueError("Model type {0} is not a valid model for simulation".format(model))

        # Copy over tree assuming units in substitutions per site
        assert isinstance(tree, Bio.Phylo.BaseTree.Tree), "invalid tree"
        self._tree = copy.deepcopy(tree)
        self._root = False
        self._internalnode = []
        self._terminalnode = []
        self._seq_array = []
        self._cached_exp_M = {}
        self._cached_cumsum = {}
        # For internal storage, branch lengths should be in model units
        # rather than codon substitutions per site. Here we adjust them from
        # For internal storage, branch lengths are in model units rather
        # codon substitutions per site to model units.
        if branchScale is None:
            branchScale = self._model.branchScale
        else:
            assert isinstance(branchScale, float) and branchScale > 0
        for node in self._tree.find_clades():
            if node != self._tree.root:
                node.branch_length /= branchScale  # adjust units
                if len(node.clades) == 0:  # terminal node
                    self._terminalnode.append(node.name)
                else:
                    self._internalnode.append(node.name)
            else:
                assert not self._root, "Tree has > 1 root. Please re-root tree"
                self._root = node
                self._internalnode.append(node.name)

        # set up sequence arrays
        for x in range(N_CODON):
            site_seq = scipy.zeros(N_CODON)
            site_seq[x] = 1
            self._seq_array.append(site_seq)
        self._seq_array = scipy.array(self._seq_array)


    def simulate(self, randomSeed=False):
        """Simulate an alignment.

        The root sequence is randomly drawn from the stationary
        state of the model. To ensure reproducible results,
        set the `scipy` random seed with the flag `randomSeed`.

        Attributes:
            `randomSeed` (`int` or `False`)
                Seed used to set the `scipy` random seed.

        Returns:
            `simulated_alignment` (`list`)
                Final codon alignment. Tuples of (seq_id, seq) for each
                terminal node in the tree.

        """
        def _evolve_branch(parent, child, alignment):
            """Generate new sequence given a parent and a child node."""
            branch_length = parent.distance(child)

            if branch_length not in self._cached_exp_M:
                self._cached_exp_M[branch_length] = self._model.M(branch_length)
            exp_M = self._cached_exp_M[branch_length]

            new_seq = []
            for r in range(self.nsites):
                parent_codon = alignment[parent.name][r]
                query = (r, parent_codon, branch_length)
                if query not in self._cached_cumsum:
                    self._cached_cumsum[query] = scipy.cumsum(self._seq_array[parent_codon].dot(exp_M[r]))
                cumsum = self._cached_cumsum[query]

                # choose from the new codon distribution for the site
                codon = scipy.argmin(cumsum < scipy.random.rand())
                new_seq.append(codon)

            alignment[child.name] = new_seq
            return alignment

        def _traverse_tree(parent, child, alignment):
            """Walk along a branch of the tree.

            Calls the `_evolve_branch` function for internal nodes.
            """
            if parent is None:  # at the root, need to generate a sequence
                root_seq = []
                for r in range(self.nsites):
                    # draw codon from stationary state
                    ss = self._model.stationarystate[r]
                    cumsum = scipy.cumsum(ss)
                    codon = scipy.argmin(cumsum / cumsum[-1] < scipy.random.rand())
                    # create sequence array. 1 indicates codon sequence
                    root_seq.append(codon)
                alignment[child.name] = scipy.array(root_seq)  # array of size `nsites`, each element array of 61
            else:  # internal branch, need to evolve along the branch
                alignment = _evolve_branch(parent, child, alignment)

            # Continue to traverse tree if child is an internal node
            if not child.is_terminal():
                for gchild in child.clades:
                    alignment = _traverse_tree(child, gchild, alignment)
            return alignment

        # beginning of `simulate` function
        if randomSeed is not False:
            scipy.random.seed(randomSeed)

        # set up alignment and begin tree traversal
        nodes = self._internalnode + self._terminalnode
        alignment = {node: [] for node in nodes}
        alignment = _traverse_tree(None, self._root, alignment)  # simulate the sequences

        # reformat the simulated alignment
        # turn the sequence arrays into codon sequnces
        # format alignment as a list of tuples
        simulated_alignment = []
        for node_name in self._terminalnode:
            seq = alignment[node_name]
            final = []
            for codon in seq:
                nt = INDEX_TO_CODON[codon]
                final.append(nt)
            final = "".join(final)
            assert (len(final) / 3) == self.nsites, "Unexpected sequence length"
            simulated_alignment.append((node_name, final))
        return simulated_alignment


# simulations with `pyvolve`
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
        createdAlignments = ['{0}_simulatedalignment.fasta'
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
