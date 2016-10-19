"""Tree likelihoods.

Nucleotides, amino acids, and codons are indexed by integers 0, 1, ...
using the indexing schemes defined in `phydmslib.constants`.
"""


import scipy
import Bio.Phylo
import phydmslib.models
from phydmslib.constants import *


class TreeLikelihood:
    """Uses alignment, model, and tree to calculate likelihoods.

    See `__init__` for how to initialize a `TreeLikelihood`.

    A `TreeLikelihood` object is designed to only work for a fixed
    tree topology. In the internal workings, this fixed topology is
    transformed into the node indexing scheme defined by `name_to_nodeindex`.

    Attributes:
        `tree` (instance of `Bio.Phylo.BaseTree.Tree` derived class)
            Phylogenetic tree. 
        `model` (instance of `phydmslib.models.Model` derived class)
            Specifies the substitution model for `nsites` codon sites.
        `alignment` (list of 2-tuples of strings, `(head, seq)`)
            Aligned protein-coding codon sequences. Headers match
            tip names in `tree`; sequences contain `nsites` codons.
        `loglik` (`float`)
            Current log likelihood.
        `siteloglik` (`numpy.ndarray` of floats, length `nsites`)
            `siteloglik[r]` is current log likelihood at site `r`.
        `nsites` (int)
            Number of codon sites.
        `nseqs` (int)
            Number of sequences in `alignment`, and tip nodes in `tree`.
        `nnodes` (int)
            Total number of nodes in `tree`, counting tip and internal.
            Node are indexed 0, 1, ..., `nnodes - 1`. This indexing is done
            such that the descendants of node `n` always have indices < `n`.
        `name_to_nodeindex` (`dict`)
            For each sequence name (these are headers in `alignmnent`, 
            tip names in `tree`), `name_to_nodeindex[name]` is the `int`
            index of that `node` in the internal indexing scheme.
        `internalnodes` (`list` of `int`, length `nnodes`)
            List of the indices of internal nodes only. Given the node 
            indexing scheme, the root node will always be the last one
            in this list.
        `rdescend` and `ldescend` (each a `list` of `int` of length `nnode`)
            For each node `n`, if `n` is an internal node then `rdescend[n]`
            is its right descendant and `ldescend[n]` is its left descendant.
            If `n` is a tip node, then `rdescend[n] == ldescend[n] == -1`.
            Note that the indexing scheme ensures that `rdescend[n] < n`
            and `ldescend[n] < n` for all `n` in `internalnodes`.
        `t` (`list` of `float`, length `nnodes - 1`)
            `t[n]` is the branch length leading node `n`. The branch leading
            to the root node (which has index `nnodes - 1`) is undefined
            and not used, which is why `t` is of length `nnodes - 1` rather
            than `nnodes`.
        `L` (`numpy.ndarray` of `float`, shape `(nnodes, nsites, N_CODON)`)
            `L[n][r][x]` is the partial conditional likelihood of codon 
            `x` at site `r` at node `n`.
    """

    def __init__(self, tree, alignment, model):
        """Initialize a `TreeLikelihood` object.

        Args:
            `tree`, `model`, `alignment`
                Attributes of same name described in class doc string.
        """
        assert isinstance(model, phydmslib.models.Model), "invalid model"
        self.model = model
        self.nsites = self.model.nsites

        assert isinstance(alignment, list) and all([isinstance(tup, tuple)
                and len(tup) == 2 for tup in alignment]), ("alignment is "
                "not list of (head, seq) 2-tuples")
        assert all([len(seq) == 3 * self.nsites for (head, seq) in alignment])
        assert set([head for (head, seq) in alignment]) == set([clade.name for
                clade in tree.get_terminals()])
        self.alignment = alignment
        self.nseqs = len(alignment)
        alignment_d = dict(self.alignment)

        assert isinstance(tree, Bio.Phylo.BaseTree.Tree), "invalid tree"
        self.tree = tree
        assert self.tree.count_terminals() == self.nseqs

        # index nodes
        self.name_to_nodeindex = {}
        self.nnodes = tree.count_terminals() + len(tree.get_nonterminals())
        self.rdescend = [-1] * self.nnodes
        self.ldescend = [-1] * self.nnodes
        self.t = [-1] * (self.nnodes - 1)
        self.internalnodes = [] 
        self.L = scipy.full((self.nnodes, self.nsites, N_CODON), -1, dtype='float')
        for (n, node) in enumerate(self.tree.find_clades(order='postorder')):
            if node != self.tree.root:
                assert n < self.nnodes - 1
                self.t[n] = node.branch_length
            if node.is_terminal():
                seq = alignment_d[node.name]
                for r in range(self.nsites):
                    codon = seq[3 * r : 3 * r + 3]
                    if codon == '---':
                        scipy.copyto(self.L[n][r], scipy.ones(N_CODON, dtype='float'))
                    elif codon in CODON_TO_INDEX:
                        scipy.copyto(self.L[n][r], scipy.zeros(N_CODON, dtype='float'))
                        self.L[n][r][CODON_TO_INDEX[codon]] = 1.0
                    else:
                        raise ValueError("Bad codon {0} in {1}".format(codon, node.name))
            else:
                assert len(node.clades) == 2, "not 2 children: {0}".format(node.name)
                self.rdescend[n] = self.name_to_nodeindex[node.clades[0]]
                self.ldescend[n] = self.name_to_nodeindex[node.clades[1]]
                self.internalnodes.append(n)
            self.name_to_nodeindex[node] = n
        assert len(self.internalnodes) == len(tree.get_nonterminals())

        # now update internal attributes related to likelihood
        self._updateInternals()

    def _updateInternals(self):
        """Update internal attributes related to likelihood.

        Should be called any time branch lengths or model parameters
        are changed.
        """
        with scipy.errstate(over='raise', under='raise', divide='raise',
                invalid='raise'):
            self._computePartialLikelihoods()
            self.siteloglik = scipy.log(scipy.sum(self.L[-1] * 
                    self.model.stationarystate, axis=1))
            self.loglik = scipy.sum(self.siteloglik)

    def _computePartialLikelihoods(self):
        """Update `L`."""
        for n in self.internalnodes:
            nright = self.rdescend[n]
            nleft = self.ldescend[n]
            tright = self.t[nright]
            tleft = self.t[nleft]
            Mright = self.model.M(tright)
            Mleft = self.model.M(tleft)
            # using this approach to broadcast matrix-vector mutiplication
            # http://stackoverflow.com/questions/26849910/numpy-matrix-multiplication-broadcast
            scipy.copyto(self.L[n], 
                    scipy.sum(Mright * self.L[nright][:, None, :], axis=2) *
                    scipy.sum(Mleft * self.L[nleft][:, None, :], axis=2))



if __name__ == '__main__':
    import doctest
    doctest.testmod()
