"""Tree likelihoods.

Nucleotides, amino acids, and codons are indexed by integers 0, 1, ...
using the indexing schemes defined in `phydmslib.constants`.
"""


import sys
import copy
import scipy
import scipy.optimize
import Bio.Phylo
import phydmslib.models
from phydmslib.numutils import *
from phydmslib.constants import *


class TreeLikelihood(object):
    """Uses alignment, model, and tree to calculate likelihoods.

    See `__init__` for how to initialize a `TreeLikelihood`.

    Most commonly, you will initiate a tree likelihood object, and
    then maximize it using the `maximizeLikelihood` method.

    A `TreeLikelihood` object is designed to only work for a fixed
    tree topology. In the internal workings, this fixed topology is
    transformed into the node indexing scheme defined by `name_to_nodeindex`.

    After initialization, most attributes should **not** be changed directly.
    They should only be altered via the `updateParams` method or
    by altering the `paramsarray` property. Both of these methods
    of updating attributes correctly update all of the dependent
    properties; setting other attributes directly does not.

    Attributes:
        `tree` (instance of `Bio.Phylo.BaseTree.Tree` derived class)
            Phylogenetic tree. Branch lengths are in units of codon
            substitutions per site for the current `model`.
        `model` (instance of `phydmslib.models.Model` derived class)
            Specifies the substitution model for `nsites` codon sites.
        `alignment` (list of 2-tuples of strings, `(head, seq)`)
            Aligned protein-coding codon sequences. Headers match
            tip names in `tree`; sequences contain `nsites` codons.
        `paramsarray` (`numpy.ndarray` of floats, 1-dimensional)
            An array of all of the free parameters. You can directly
            assign to `paramsarray`, and the correct parameters will
            be internally updated vi `updateParams`. 
            `paramsarray` is actually a property rather than an attribute.
        `paramsarraybounds` (`tuple` of 2-tuples)
            `paramsarraybounds[i]` is the 2-tuple `(minbound, maxbound)`
            for parameter `paramsarray[i]`. Set `minbound` or `maxbound`
            are `None` if there is no lower or upper bound.
        `loglik` (`float`)
            Current log likelihood.
        `siteloglik` (`numpy.ndarray` of floats, length `nsites`)
            `siteloglik[r]` is current log likelihood at site `r`.
        `dloglik` (`dict`)
            For each `param` in `model.freeparams`, `dloglik[param]`
            is the derivative of `loglik` with respect to `param`.
        `dloglikarray` (`numpy.ndarray` of floats, 1-dimensional)
            `dloglikarray[i]` is the derivative of `loglik` with respect
            to the parameter represented in `paramsarray[i]`. This is
            the same information as in `dloglik`, but in a different
            representation.
        `dsiteloglik` (`dict`)
            For each `param` in `model.freeparams`, `dsiteloglik[param][r]`
            is the derivative of `siteloglik[r]` with respect to to `param`.
        `nsites` (int)
            Number of codon sites.
        `nseqs` (int)
            Number of sequences in `alignment`, and tip nodes in `tree`.
        `nnodes` (int)
            Total number of nodes in `tree`, counting tip and internal.
            Node are indexed 0, 1, ..., `nnodes - 1`. This indexing is done
            such that the descendants of node `n` always have indices < `n`,
            and so that the tips are the first `ntips` indices and the 
            internals are the last `ninternal` indices.
        `ninternal` (int)
            Total number of internal nodes in `tree`.
        `ntips` (int)
            Total number of tip nodes in `tree`.
        `tips` (`numpy.ndarray` of `int`, shape `(ntips, nsites)`)
            `tips[n][r]` is the index of the codon at site `r`
            for tip node `n` (0 <= `n` < `ntips`). If `tips[n][r]` is 
            a gap, then `tips[n][r]` is 0 and `r` is in `gaps[n]`.
        `gaps` (`numpy.array`, length `ntips`, entries `numpy.array` of `int`)
            `gaps[n]` is an array containing all sites where the codon is
            a gap for tip node `n` (0 <= `n` < `ntips`).
        `name_to_nodeindex` (`dict`)
            For each sequence name (these are headers in `alignment`, 
            tip names in `tree`), `name_to_nodeindex[name]` is the `int`
            index of that `node` in the internal indexing scheme.
        `rdescend` and `ldescend` (lists of `int` of length `ninternal`)
            For each internal node `n`, `rdescend[n - ntips]`
            is its right descendant and `ldescend[n - ntips]` is its left 
            descendant.
        `t` (`list` of `float`, length `nnodes - 1`)
            `t[n]` is the branch length leading node `n`. The branch leading
            to the root node (which has index `nnodes - 1`) is undefined
            and not used, which is why `t` is of length `nnodes - 1` rather
            than `nnodes`.
        `L` (`numpy.ndarray` of `float`, shape `(ninternal, nsites, N_CODON)`)
            `L[n - ntips][r][x]` is the partial conditional likelihood of
            codon `x` at site `r` at internal node `n`. Note that these
            must be corrected by adding `underflowlogscale`.
        `dL` (`dict` keyed by strings, values `numpy.ndarray` of `float`)
            For each free model parameter `param` in `model.freeparam`, 
            `dL[param]` is derivative of `L` with respect to `param`.
            If `param` is float, `dL[param][n - ntips][r][x]` is derivative
            of `L[n - ntips][r][x]` with respect to `param`. If `param` is 
            array, `dL[param][n - ntips][i][r][x]` is derivative of 
            `L[n][r][x]` with respect to `param[i]`.
            with respect to `param[i]`.
        `underflowfreq` (`int` >= 1)
            The frequency with which we rescale likelihoods to avoid
            numerical underflow
        `underflowlogscale` (`numpy.ndarray` of `float`, shape `(nsites,)`)
            Corrections that must be added to `L` at the root node to
            get actual likelihoods when underflow correction is being
            performed.
    """

    def __init__(self, tree, alignment, model, underflowfreq=5):
        """Initialize a `TreeLikelihood` object.

        Args:
            `tree`, `model`, `alignment`, `underflowfreq`
                Attributes of same name described in class doc string.
                Note that we make copies of both `tree` and `model`
                so the calling objects are not modified during 
                optimization.
        """
        assert isinstance(underflowfreq, int) and underflowfreq >= 1
        self.underflowfreq = underflowfreq

        self._checkModel(model)
        self.model = copy.deepcopy(model)
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

        # For internal storage, branch lengths are in model units rather
        # than codon substitutions per site. So here adjust them from
        # codon substitutions per site to model units.
        assert isinstance(tree, Bio.Phylo.BaseTree.Tree), "invalid tree"
        self._tree = copy.deepcopy(tree)
        assert self._tree.count_terminals() == self.nseqs
        branchScale = self.model.branchScale
        for node in self._tree.find_clades():
            if node != self._tree.root:
                node.branch_length /= branchScale

        # index nodes
        self.name_to_nodeindex = {}
        self.ninternal = len(self._tree.get_nonterminals())
        self.ntips = self._tree.count_terminals()
        self.nnodes = self.ntips + self.ninternal
        self.rdescend = [-1] * self.ninternal
        self.ldescend = [-1] * self.ninternal
        self.t = [-1] * (self.nnodes - 1)
        self.L = scipy.full((self.ninternal, self.nsites, N_CODON), -1, 
                dtype='float')
        self.underflowlogscale = scipy.zeros(self.nsites, dtype='float')
        self.dL = {}
        for param in self.model.freeparams:
            paramvalue = getattr(self.model, param)
            if isinstance(paramvalue, float):
                self.dL[param] = scipy.full((self.ninternal, self.nsites,
                        N_CODON), -1, dtype='float')
            elif isinstance(paramvalue, scipy.ndarray) and (paramvalue.ndim
                    == 1):
                self.dL[param] = scipy.full((self.ninternal, len(paramvalue),
                        self.nsites, N_CODON), -1, dtype='float')
            else:
                raise ValueError("Cannot handle param: {0}, {1}".format(
                        param, paramvalue))
        tipnodes = []
        internalnodes = []
        self.tips = scipy.zeros((self.ntips, self.nsites), dtype='int')
        self.gaps = []
        for node in self._tree.find_clades(order='postorder'):
            if node.is_terminal():
                tipnodes.append(node)
            else:
                internalnodes.append(node)
        for (n, node) in enumerate(tipnodes + internalnodes):
            if node != self._tree.root:
                assert n < self.nnodes - 1
                self.t[n] = node.branch_length
            if node.is_terminal():
                assert n < self.ntips
                seq = alignment_d[node.name]
                nodegaps = []
                for r in range(self.nsites):
                    codon = seq[3 * r : 3 * r + 3]
                    if codon == '---':
                        nodegaps.append(r)
                    elif codon in CODON_TO_INDEX:
                        self.tips[n][r] = CODON_TO_INDEX[codon]
                    else:
                        raise ValueError("Bad codon {0} in {1}".format(codon, 
                                node.name))
                self.gaps.append(scipy.array(nodegaps, dtype='int'))
            else:
                assert n >= self.ntips, "n = {0}, ntips = {1}".format(
                        n, self.ntips)
                assert len(node.clades) == 2, ("not 2 children: {0} has {1}\n"
                        "Is this the root node? -- {2}\n"
                        "Try `tree.root_at_midpoint()` before passing "
                        "to `TreeLikelihood`.").format(
                        node.name, len(node.clades), node == self._tree.root)
                ni = n - self.ntips
                self.rdescend[ni] = self.name_to_nodeindex[node.clades[0]]
                self.ldescend[ni] = self.name_to_nodeindex[node.clades[1]]
            self.name_to_nodeindex[node] = n
        self.gaps = scipy.array(self.gaps)
        assert len(self.gaps) == self.ntips

        # _index_to_param defines internal mapping of
        # `paramsarray` indices and parameters
        self._paramsarray = None # will be set by `paramsarray` property as needed
        self._index_to_param = {} # value is parameter associated with index
        i = 0
        for param in self.model.freeparams:
            paramvalue = getattr(self.model, param)
            if isinstance(paramvalue, float):
                self._index_to_param[i] = param
                i += 1
            else:
                for (pindex, pvalue) in enumerate(paramvalue):
                    self._index_to_param[i] = (param, pindex)
                    i += 1
        self._param_to_index = dict([(y, x) for (x, y) in 
                self._index_to_param.items()])
        assert i == len(self._index_to_param) 

        # now update internal attributes related to likelihood
        self._updateInternals()

    def _checkModel(self, model):
        """Makes sure `model` is appropriate."""
        assert isinstance(model, phydmslib.models.Model), "invalid model"
        assert not isinstance(model, phydmslib.models.DistributionModel), (
                "Cannot use a simple `Model` for `TreeLikelihood`. Use "
                "`TreeLikelihoodDistributionModel` instead if you have "
                "a `DistributionModel`.")

    def maximizeLikelihood(self, approx_grad=False):
        """Maximize the log likelihood.

        Uses the L-BFGS-B method implemented in `scipy.optimize`.

        There is no return variable, but after call object attributes
        will correspond to maximimum likelihood values.

        Args:
            `approx_grad` (bool)
                If `True`, then we numerically approximate the gradient
                rather than using the analytical values.

        Returns:
            A `scipy.optimize.OptimizeResult` with result of maximization.
        """
        # Some useful notes on optimization:
        # http://www.scipy-lectures.org/advanced/mathematical_optimization/

        def func(x):
            """Returns negative log likelihood when `x` is parameter array."""
            self.paramsarray = x
            return -self.loglik

        def dfunc(x):
            """Returns negative gradient log likelihood."""
            self.paramsarray = x
            return -self.dloglikarray

        if approx_grad:
            dfunc = False

        assert len(self.paramsarray) > 0, "There are no parameters to optimize"
        result = scipy.optimize.minimize(func, self.paramsarray,
                method='L-BFGS-B', jac=dfunc, bounds=self.paramsarraybounds)
        assert result.success, ("Optimization failed with this message:\n{0}\n"
                "Current mu = {1}\nCurrent model params: {2}").format(
                result.message, self.model.mu, ', '.join(['{0} = {1}'.format(
                p, pvalue) for (p, pvalue) in self.model.paramsReport.items()])) 

        return result

    @property
    def tree(self):
        """Tree with branch lengths in codon substitutions per site.

        The tree is a `Bio.Phylo.BaseTree.Tree` object.

        This is the current tree after whatever optimizations have
        been performed so far.
        """
        bs = self.model.branchScale
        for node in self._tree.find_clades():
            if node != self._tree.root:
                node.branch_length = self.t[self.name_to_nodeindex[node]] * bs
        return self._tree

    @property
    def paramsarraybounds(self):
        """Bounds for parameters in `paramsarray`."""
        bounds = []
        for (i, param) in self._index_to_param.items():
            if isinstance(param, str):
                bounds.append(self.model.PARAMLIMITS[param])
            elif isinstance(param, tuple):
                bounds.append(self.model.PARAMLIMITS[param[0]])
            else:
                raise ValueError("Invalid param type")
        bounds = [(tup[0] + ALMOST_ZERO, tup[1] - ALMOST_ZERO) for tup in bounds]
        assert len(bounds) == len(self._index_to_param)
        return tuple(bounds)

    @property
    def paramsarray(self):
        """All free model parameters as 1-dimensional `numpy.ndarray`.
        
        You are allowed to update model parameters by direct
        assignment of this property."""
        # Always return copy of `_paramsarray` because setter checks if changed
        if self._paramsarray is not None:
            return self._paramsarray.copy() 
        nparams = len(self._index_to_param)
        self._paramsarray = scipy.ndarray(shape=(nparams,), dtype='float')
        for (i, param) in self._index_to_param.items():
            if isinstance(param, str):
                self._paramsarray[i] = getattr(self.model, param)
            elif isinstance(param, tuple):
                self._paramsarray[i] = getattr(self.model, param[0])[param[1]]
            else:
                raise ValueError("Invalid param type")
        return self._paramsarray.copy()

    @paramsarray.setter
    def paramsarray(self, value):
        """Set new `paramsarray` and update via `updateParams`."""
        nparams = len(self._index_to_param)
        assert (isinstance(value, scipy.ndarray) and value.ndim == 1), (
                "paramsarray must be 1-dim ndarray")
        assert len(value) == nparams, ("Assigning paramsarray to ndarray "
                "of the wrong length.")
        if (self._paramsarray is not None) and all(value == self._paramsarray):
            return # do not need to do anything if value has not changed
        # build `newvalues` to pass to `updateParams`
        newvalues = {}
        vectorized_params = {}
        for (i, param) in self._index_to_param.items():
            if isinstance(param, str):
                newvalues[param] = float(value[i])
            elif isinstance(param, tuple):
                (iparam, iparamindex) = param
                if iparam in vectorized_params:
                    assert iparamindex not in vectorized_params[iparam]
                    vectorized_params[iparam][iparamindex] = float(value[i])
                else:
                    vectorized_params[iparam] = {iparamindex:float(value[i])}
            else:
                raise ValueError("Invalid param type")
        for (param, paramd) in vectorized_params.items():
            assert set(paramd.keys()) == set(range(len(paramd)))
            newvalues[param] = scipy.array([paramd[i] for i in range(len(paramd))],
                    dtype='float')
        self.updateParams(newvalues)
        self._paramsarray = self.paramsarray

    @property
    def dloglikarray(self):
        """Derivative of `loglik` with respect to `paramsarray`."""
        nparams = len(self._index_to_param)
        dloglikarray = scipy.ndarray(shape=(nparams,), dtype='float')
        for (i, param) in self._index_to_param.items():
            if isinstance(param, str):
                dloglikarray[i] = self.dloglik[param]
            elif isinstance(param, tuple):
                dloglikarray[i] = self.dloglik[param[0]][param[1]]
        return dloglikarray

    def updateParams(self, newvalues):
        """Update parameters and re-compute likelihoods.

        This method is the **only** acceptable way to update model
        or tree parameters. The likelihood is re-computed as needed
        by this method.

        Args:
            `newvalues` (dict)
                A dictionary keyed by param name and with value as new
                value to set. Each parameter name must either be
                a string in `model.freeparams` or a string that represents
                a valid branch length.
        """
        modelparams = {}
        otherparams = {}
        for (param, value) in newvalues.items():
            if param in self.model.freeparams:
                modelparams[param] = value
            else:
                otherparams[param] = value
        if modelparams:
            self.model.updateParams(modelparams)
        if otherparams:
            raise RuntimeError("Cannot handle non-model params: {0}".format(
                    otherparams))
        if newvalues:
            self._updateInternals()
            self._paramsarray = None

    def _updateInternals(self):
        """Update internal attributes related to likelihood.

        Should be called any time branch lengths or model parameters
        are changed.
        """
        with scipy.errstate(over='raise', under='raise', divide='raise',
                invalid='raise'):
            self.underflowlogscale.fill(0.0)
            self._computePartialLikelihoods()
            sitelik = scipy.sum(self.L[-1] * self.model.stationarystate, axis=1)
            self.siteloglik = scipy.log(sitelik) + self.underflowlogscale
            self.loglik = scipy.sum(self.siteloglik)
            self.dsiteloglik = {}
            self.dloglik = {}
            for param in self.model.freeparams:
                self.dsiteloglik[param] = scipy.sum(self.model.dstationarystate(param)
                        * self.L[-1] + self.dL[param][-1] 
                        * self.model.stationarystate, axis=-1) / sitelik
                self.dloglik[param] = scipy.sum(self.dsiteloglik[param], axis=-1)

    def _computePartialLikelihoods(self):
        """Update `L`. and `dL`."""
        for n in range(self.ntips, self.nnodes):
            ni = n - self.ntips # internal node number
            nright = self.rdescend[ni]
            nleft = self.ldescend[ni]
            nrighti = nright - self.ntips # internal node number
            nlefti = nleft - self.ntips # internal node number
            if nright < self.ntips:
                istipr = True
            else:
                istipr = False
            if nleft < self.ntips:
                istipl = True
            else:
                istipl = False
            tright = self.t[nright]
            tleft = self.t[nleft]
            if istipr:
                Mright = MLright = self.model.M(tright, self.tips[nright], 
                        self.gaps[nright])
            else:
                Mright = self.model.M(tright)
                MLright = broadcastMatrixVectorMultiply(Mright, 
                        self.L[nrighti])
            if istipl:
                Mleft = MLleft = self.model.M(tleft, self.tips[nleft], 
                        self.gaps[nleft])
            else:
                Mleft = self.model.M(tleft)
                MLleft = broadcastMatrixVectorMultiply(Mleft, self.L[nlefti])
            scipy.copyto(self.L[ni], MLright * MLleft)
            if ni > 0 and ni % self.underflowfreq == 0:
                scale = scipy.amax(self.L[ni], axis=1)
                self.L[ni] /= scale[:, scipy.newaxis]
                self.underflowlogscale += scipy.log(scale)
            for param in self.model.freeparams:
                paramvalue = getattr(self.model, param)
                if isinstance(paramvalue, float):
                    indices = [()] # no sub-indexing needed
                elif (isinstance(paramvalue, scipy.ndarray) and 
                        paramvalue.ndim == 1 and paramvalue.shape[0] > 1):
                    # need to sub-index calculations for each param element
                    indices = [(j,) for j in range(len(paramvalue))]
                else:
                    raise RuntimeError("invalid param: {0}".format(param))
                if istipr:
                    dMright = self.model.dM(tright, param, Mright, 
                            self.tips[nright], self.gaps[nright])
                else:
                    dMright = self.model.dM(tright, param, Mright)
                if istipl:
                    dMleft = self.model.dM(tleft, param, Mleft, 
                            self.tips[nleft], self.gaps[nleft])
                else:
                    dMleft = self.model.dM(tleft, param, Mleft)
                for j in indices:
                    if istipr:
                        dMLright = dMright[j]
                        MdLright = 0
                    else:
                        dMLright = broadcastMatrixVectorMultiply(
                                dMright[j], self.L[nrighti])
                        MdLright = broadcastMatrixVectorMultiply(Mright,
                                self.dL[param][nrighti][j])
                    if istipl:
                        dMLleft = dMleft[j]
                        MdLleft = 0
                    else:
                        dMLleft = broadcastMatrixVectorMultiply(
                                dMleft[j], self.L[nlefti])
                        MdLleft = broadcastMatrixVectorMultiply(Mleft,
                                self.dL[param][nlefti][j])
                    scipy.copyto(self.dL[param][ni][j], (dMLright + MdLright)
                            * MLleft + MLright * (dMLleft + MdLleft))
                    if ni > 0 and ni % self.underflowfreq == 0:
                        self.dL[param][ni][j] /= scale[:, scipy.newaxis] 


class TreeLikelihoodDistribution(TreeLikelihood):
    """Like `TreeLikelihood` but for `DistributionModel`.
   
    Differs from `TreeLikelihood` in that `model` is a `DistributionModel`
    rather than just a simple `Model`. Otherwise inherits all attributes
    and functions from `TreeLikelihood` with the following differences:

    Attributes that differ from `TreeLikelihood` base model:
        `model` (instance of `phydmslib.models.DistributionModel`)
        `ncats` (`int`)
            Number of categories for distributed parameter of `model`.
            Determined from `model`.
        `siteloglik` (`numpy.ndarray` of float, shape `(ncats, nsites)`)
            `siteloglik[k][r]` is log likelihood at site `r` for
            category `k`.
        `dsiteloglik` (`dict`)
            For each free parameter in `model.freeparams` that is 
            **not** also in `model.distributionparams`, plus for
            `model.distributedparam`, `dsiteloglik[param][k][r]` is
            the derivative of `siteloglik[k][r]` with respect to
            `param`.
        `L` (`numpy.ndarray`, shape `(ninternal, ncats, nsites, N_CODON)`)
            `L[n - ntips][k][r][x]` is partial conditional likelihood
            of `x` at `r` at internal node `n` for category `k`.
            Note that these must be corrected by adding 
            `underflowlogscale`.
        `dL` (`dict` keyed by strings)
            For each free parameter that keys `dsiteloglik`,
            `dL[param][n - ntips][k]` is the derivative of
            `L[n - ntips][k]` with respect to `param`.
    """

    def __init__(self, tree, alignment, model, underflowfreq=5):
        """Initialize a `TreeLikelihoodDistribution` object.

        See docs for `TreeLikelihood.__init__`."""

        super(TreeLikelihoodDistribution, self).__init__(
                tree, alignment, model, underflowfreq=underflowfreq)

        self.ncats = self.model.ncats

    def _resize_attributes(self):
        """Re-size attributes for `TreeLikelihoodDistribution`.

        Some of the attributes for a `TreeLikelihoodDistribution`
        have a different size than for a `TreeLikelihood` due to
        need to store data for the `ncats` categories. This method
        checks the size of all such attributes, and if necessary
        re-sizes them.

        This is needed because the `super` call to `__init__`
        will leave sizes for `TreeLikelihood` rather than
        `TreeLikelihoodDistribution`.
        """
        Lshape = (self.ninternal, self.ncats, self.nsites, N_CODON)
        if self.L.shape != Lshape:
            self.L = scipy.full(Lshape, -1, dtype='float')
        for (param, paramdL) in list(self.dL.item()):
            if paramdL.ndim == 3:
                self.dL[param] = scipy.full(Lshape, -1, dtype='float')
            elif paramdL.ndim == 4:
                newshape = (self.ninternal, self.ncats, 
                        paramdL.shape[1], self.nsites, N_CODON)
                self.dL[param] = scipy.full(newshape, -1, 
                        dtype='float')
            else:
                raise RuntimeError("Invalid shape: {0}, {1}".format(
                        param, paramdL.shape)

    def _checkModel(self, model):
        """Makes sure `model` is appropriate."""
        assert isinstance(model, phydmslib.models.DistributionModel), (
                "TreeLikelihoodDistribution requires DistributionModel")

    def _updateInternals(self):
        """Update internal attributes related to likelihood.

        Should be called anytime branch lengths or model parameters
        are changed.
        """
        self._resize_attributes()
        raise RuntimeError('not yet implemented')

    def _computePartialLikelihoods(self):
        """Update `L` and `dL`."""
        raise RuntimeError('not yet implemented')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
