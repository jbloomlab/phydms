"""Tree likelihoods.

Nucleotides, amino acids, and codons are indexed by integers 0, 1, ...
using the indexing schemes defined in `phydmslib.constants`.
"""


import sys
import math
import copy
import warnings
warnings.simplefilter('always')
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
            This can either be a simple `Model` or a `DistributionModel`.
            Note that if a non-zero prior is defined by `model`, this
            is included in `loglik` and `dloglik`.
        `alignment` (list of 2-tuples of strings, `(head, seq)`)
            Aligned protein-coding codon sequences. Headers match
            tip names in `tree`; sequences contain `nsites` codons.
        `paramsarray` (`numpy.ndarray` of floats, 1-dimensional)
            An array of the model free parameters. You can directly
            assign to `paramsarray`, and the correct parameters will
            be internally updated vi `updateParams`.
            `paramsarray` is actually a property rather than an attribute.
        `paramsarraybounds` (`tuple` of 2-tuples)
            `paramsarraybounds[i]` is the 2-tuple `(minbound, maxbound)`
            for parameter `paramsarray[i]`. Set `minbound` or `maxbound`
            are `None` if there is no lower or upper bound.
        `loglik` (`float`)
            Current log likelihood. Also includes `model.logprior`
            if this is non-zero.
        `siteloglik` (`numpy.ndarray` of floats, length `nsites`)
            `siteloglik[r]` is current log likelihood at site `r`.
        `dparamscurrent` (`bool`)
            Do we keep the derivatives of the likelihood with
            respect to model parameters current? Doing so
            involves computational costs, so only set to `True`
            if using these derivatives. Currently, only one of
            `dparamscurrent` and `dtcurrent` can be `True` simultaneously,
            as we optimize parameters and branch lengths separately.
        `dtcurrent` (`bool`)
            Do we keep the derivatives of the likelihood with
            respect to branch lengths current? Doing so
            involves computational costs, so only set to `True`
            if using these derivatives.
        `dloglik` (`dict`)
            For each `param` in `model.freeparams`, `dloglik[param]`
            is the derivative of `loglik` with respect to `param`.
            Only accessible if `dparamscurrent` is `True`. Also
            includes `model.dlogprior` if this is non-zero.
        `dloglikarray` (`numpy.ndarray` of floats, 1-dimensional)
            `dloglikarray[i]` is the derivative of `loglik` with respect
            to the parameter represented in `paramsarray[i]`. This is
            the same information as in `dloglik`, but in a different
            representation. Only accessible if `dparamscurrent` is `True`.
        `dloglik_dt` (`numpy.ndarray` of floats, shape `(nnodes - 1,)`.
            `dloglik_dt` is the derivative of `loglik` with respect
            to the branch lengths `t`.
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
        `descendants` (`list` of `list`)
            For each node `n`, `descendants[n]` is a list of all nodes
            that are direct or indirect descendants of `n`.
        `t` (`numpy.ndarray` of `float`, shape `(nnodes - 1,)`)
            `t[n]` is the branch length leading node `n`. The branch leading
            to the root node (which has index `nnodes - 1`) is undefined
            and not used, which is why `t` is of length `nnodes - 1` rather
            than `nnodes`.
        `L` (`dict`)
            `L[n]` is a `numpy.ndarray` of `float` of shape `(nsites, N_CODON)`.
            for each node `n`. `L[n][r][x] is the partial condition likelihood
            of codon `x` at site `r` at node `n`. Note that these
            must be corrected by adding `underflowlogscale`.
            If `model` is a `DistributionModel`, then `L[n]` is instead
            of shape `(model.ncats, nsites, N_CODON)` with the
            first index ranging over the model categories.
            Entries in `L` are only defined and saved for nodes as needed
            in order to save memory.
        `dL` (`dict` keyed by strings)
            For each free model parameter `param` and node `n`,
            `dL[param][n]` is derivative of `L[n]` with respect to `param`.
            Only guaranteed to be valid if `dparamscurrent` is `True`,
            and entries are only defined and saved for nodes as needed
            in order to save memory.
        `dL_dt` (`dict` of `dict`)
            `dL_dt[n][n2]` is the derivative of `L[n2]` with respect to
            branch length `t[n]` where `0 <= n < nnodes - 1`.
        `underflowfreq` (`int` >= 1)
            The frequency with which we rescale likelihoods to avoid
            numerical underflow
        `underflowlogscale` (`numpy.ndarray` of `float`, shape `(nsites,)`)
            Corrections that must be added to `L` at the root node to
            get actual likelihoods when underflow correction is being
            performed.
    """

    def __init__(self, tree, alignment, model, underflowfreq=3,
            dparamscurrent=True, dtcurrent=False, branchScale=None):
        """Initialize a `TreeLikelihood` object.

        Args:
            `tree`, `alignment`, `model`, `underflowfreq`
                Attributes of same name described in class doc string.
                Note that we make copies of `tree`, `model`, and
                `alignment`' so the calling objects are not modified
                during optimization.
            `dparamscurrent`, `dtcurrent`
                Current values of these flags indicating which
                derivatives to compute.
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
        assert isinstance(underflowfreq, int) and underflowfreq >= 1
        self.underflowfreq = underflowfreq

        assert isinstance(dparamscurrent, bool)
        self._dparamscurrent = dparamscurrent
        assert isinstance(dtcurrent, bool)
        self._dtcurrent = dtcurrent
        assert not (self.dparamscurrent and self.dtcurrent), (
                "Only one of dparamscurrent or dtcurrent can be True")

        if isinstance(model, phydmslib.models.DistributionModel):
            self._distributionmodel = True
        elif isinstance(model, phydmslib.models.Model):
            self._distributionmodel = False
        else:
            raise ValueError("model is not a valid model")
        self.model = copy.deepcopy(model)
        self.nsites = self.model.nsites

        assert isinstance(alignment, list) and all([isinstance(tup, tuple)
                and len(tup) == 2 for tup in alignment]), ("alignment is "
                "not list of (head, seq) 2-tuples")
        assert all([len(seq) == 3 * self.nsites for (head, seq) in alignment])
        assert set([head for (head, seq) in alignment]) == set([clade.name for
                clade in tree.get_terminals()])
        self.alignment = copy.deepcopy(alignment)
        self.nseqs = len(alignment)
        alignment_d = dict(self.alignment)

        # For internal storage, branch lengths are in model units rather
        # than codon substitutions per site. So here adjust them from
        # codon substitutions per site to model units.
        assert isinstance(tree, Bio.Phylo.BaseTree.Tree), "invalid tree"
        self._tree = copy.deepcopy(tree)
        assert self._tree.count_terminals() == self.nseqs
        if branchScale == None:
            branchScale = self.model.branchScale
        else:
            assert isinstance(branchScale, float) and branchScale > 0
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
        self.underflowlogscale = scipy.zeros(self.nsites, dtype='float')
        if self._distributionmodel:
            self._Lshape = (self.model.ncats, self.nsites, N_CODON)
        else:
            self._Lshape = (self.nsites, N_CODON)
        self.dL_dt = dict([(n, {}) for n in range(self.nnodes - 1)])
        self.L = {}
        self.dL = {}
        self._dLshape = {}
        for param in self._paramlist_PartialLikelihoods:
            self.dL[param] = {}
            if self._distributionmodel and param == self.model.distributedparam:
                self._dLshape[param] = self._Lshape
            else:
                paramvalue = getattr(self.model, param)
                if isinstance(paramvalue, float):
                    self._dLshape[param] = self._Lshape
                elif isinstance(paramvalue, scipy.ndarray) and (paramvalue.ndim
                        == 1):
                    if self._distributionmodel:
                        self._dLshape[param] = (self.model.ncats,
                                len(paramvalue), self.nsites, N_CODON)
                    else:
                        self._dLshape[param] = (len(paramvalue),
                                self.nsites, N_CODON)
                else:
                    raise ValueError("Cannot handle param: {0}, {1}".format(
                            param, paramvalue))
        tipnodes = []
        internalnodes = []
        self.tips = scipy.zeros((self.ntips, self.nsites), dtype='int')
        self.gaps = []
        self.descendants = []
        self._t = scipy.full((self.nnodes - 1,), -1, dtype='float')
        for node in self._tree.find_clades(order='postorder'):
            if node.is_terminal():
                tipnodes.append(node)
            else:
                internalnodes.append(node)
        for (n, node) in enumerate(tipnodes + internalnodes):
            if node != self._tree.root:
                assert n < self.nnodes - 1
                self._t[n] = node.branch_length
            if node.is_terminal():
                assert n < self.ntips
                self.descendants.append([])
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
                self.descendants.append([self.name_to_nodeindex[nx] for nx
                        in node.find_clades() if nx != node])
            self.name_to_nodeindex[node] = n
        self.gaps = scipy.array(self.gaps)
        assert len(self.gaps) == self.ntips

        # _index_to_param defines internal mapping of
        # `paramsarray` indices and parameters
        self._paramsarray = None # set by `paramsarray` property as needed
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

    def maximizeLikelihood(self, optimize_brlen=False,
            approx_grad=False, logliktol=1.0e-2, nparamsretry=1,
            printfunc=None):
        """Maximize the log likelihood.

        Maximizes log likelihood with respect to model parameters
        and potentially branch lengths depending on `optimize_brlen`.
        If optimizing the branch lengths, iterates between optimizing
        the model parameters and branch lengths.

        Uses the L-BFGS-B method implemented in `scipy.optimize`.

        There is no return variable, but after call object attributes
        will correspond to maximimum likelihood values.

        Args:
            `optimize_brlen` (bool)
                Do we optimize branch lengths?
            `approx_grad` (bool)
                If `True`, then we numerically approximate the gradient
                rather than using the analytical values.
            `logliktol` (float)
                When using `optimize_brlen`, keep iterating between
                optimization of parameters and branch lengths until
                change in log likelihood is less than `logliktol`.
            `nparamsretry` (int >= 0)
                Number of times to retry parameter optimization from
                different initial values if it fails the first time.
            `printfunc` (`None` or a function)
                If not `None`, then we print using `printfunc` the
                detailed results of the optimization at each step.
                For instance, `printfunc` might be `sys.stderr.write`
                or `logger.info`.

        Returns:
            A string giving a summary of the maximization.
        """
        # Some useful notes on optimization:
        # http://www.scipy-lectures.org/advanced/mathematical_optimization/

        assert len(self.paramsarray) > 0, "No parameters to optimize"
        assert nparamsretry >= 0
        assert logliktol > 0

        def paramsfunc(x):
            """Negative log likelihood when `x` is params."""
            self.paramsarray = x
            return -self.loglik

        def paramsdfunc(x):
            """Negative gradient log likelihood with respect to params."""
            self.paramsarray = x
            return -self.dloglikarray

        def tfunc(x):
            """Negative log likelihood when `x` is branch lengths."""
            self.t = x
            return -self.loglik

        def tdfunc(x):
            """Negative gradient loglik with respect to branch lengths."""
            self.t = x
            return -self.dloglik_dt

        if approx_grad:
            paramsdfunc = False
            tdfunc = False
            self.dtcurrent = False
            self.dparamscurrent = False

        def _printResult(opttype, result, i, old, new):
            """Print summary of optimization result."""
            if printfunc is not None:
                printfunc('Step {0}, optimized {1}.\n'
                          'Likelihood went from {2} to {3}.\n'
                          'Max magnitude in Jacobian is {4}.\n'
                          'Full optimization result:\n{5}\n'.format(
                          i, opttype, old, new,
                          scipy.absolute(result.jac).max(), result))

        oldloglik = self.loglik
        converged = False
        firstbrlenpass = True
        options = {'ftol':1.0e-7} # optimization options
        summary = []
        i = 1
        while not converged:
            if (not self.dparamscurrent) and (not approx_grad):
                self.dtcurrent = False
                self.dparamscurrent = True
            nparamstry = 0
            origparamsarray = self.paramsarray.copy()
            paramsconverged = False
            while not paramsconverged:
                result = scipy.optimize.minimize(paramsfunc, self.paramsarray,
                        method='L-BFGS-B', jac=paramsdfunc,
                        bounds=self.paramsarraybounds, options=options)
                _printResult('params', result, i, oldloglik, self.loglik)
                msg = ('Step {0}: optimized parameters, loglik went from '
                        '{1:.2f} to {2:.2f} ({3} iterations, {4} function '
                        'evals)'.format(i, oldloglik, self.loglik, result.nit,
                        result.nfev))
                summary.append(msg)
                if result.success and (not (oldloglik - self.loglik > logliktol)):
                    paramsconverged = True
                    jacmax = scipy.absolute(result.jac).max()
                    if (jacmax > 1000) and not (firstbrlenpass and optimize_brlen):
                        warnings.warn("Optimizer reports convergence, "
                                "but max element in Jacobian is {0}\n"
                                "Summary of optimization:\n{1}".format(
                                jacmax, summary))
                else:
                    if not result.success:
                        resultmessage = result.message
                    else:
                        resultmessage = ('loglik increased in param optimization '
                                'from {0} to {1}'.format(oldloglik, self.loglik))
                    nparamstry += 1
                    failmsg = ("Optimization failure {0}\n{1}\n{2}".format(
                            nparamstry, resultmessage, '\n'.join(summary)))
                    if nparamstry > nparamsretry:
                        raise RuntimeError(failmsg)
                    else:
                        warnings.warn(failmsg + '\n\n' +
                                "Re-trying with different initial params.")
                        scipy.random.seed(nparamstry)
                        # seed at geometric mean of original value, max
                        # bound, min bound, and random number between max and min
                        minarray = scipy.array([self.paramsarraybounds[j][0] for
                                j in range(len(self.paramsarray))])
                        maxarray = scipy.array([self.paramsarraybounds[j][1] for
                                j in range(len(self.paramsarray))])
                        randarray = scipy.random.uniform(minarray, maxarray)
                        newarray = (minarray * maxarray * randarray *
                                origparamsarray)**(1 / 4.) # geometric mean
                        assert newarray.shape == self.paramsarray.shape
                        assert (newarray > minarray).all()
                        assert (newarray < maxarray).all()
                        self.paramsarray = newarray
            i += 1
            assert oldloglik - self.loglik <= logliktol
            if (self.loglik - oldloglik >= logliktol) or firstbrlenpass:
                firstbrlenpass = False
                oldloglik = self.loglik
                if optimize_brlen:
                    if not approx_grad:
                        self.dparamscurrent = False
                        self.dtcurrent = True
                    result = scipy.optimize.minimize(tfunc, self.t,
                            method='L-BFGS-B', jac=tdfunc, options=options,
                            bounds=[(ALMOST_ZERO, None)] * len(self.t))
                    _printResult('branches', result, i, oldloglik, self.loglik)
                    summary.append('Step {0}: optimized branches, loglik '
                            'went from {1:.2f} to {2:.2f} ({3} iterations, '
                            '{4} function evals)'.format(i, oldloglik,
                            self.loglik, result.nit, result.nfev))
                    i += 1
                    assert result.success, ("Optimization failed\n{0}"
                            "\n{1}\n{2}".format(result.message, self.t,
                            '\n'.join(summary)))
                    if oldloglik - self.loglik > logliktol:
                        raise RuntimeError("loglik increased during t "
                                "optimization: {0} to {1}".format(
                                oldloglik, self.loglik))
                    elif self.loglik - oldloglik >= logliktol:
                        oldloglik = self.loglik
                    else:
                        converged = True
                else:
                    converged = True
            else:
                converged = True

        return '\n'.join(summary)

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
        # Return copy of `_paramsarray` because setter checks if changed
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
    def dparamscurrent(self):
        """Are derivatives with respect to model parameters current?"""
        return self._dparamscurrent

    @dparamscurrent.setter
    def dparamscurrent(self, value):
        """Set value of `dparamscurrent`, update derivatives if needed."""
        assert isinstance(value, bool)
        if value and self.dtcurrent:
            raise RuntimeError("Can't set both dparamscurrent and dtcurrent True")
        if value != self.dparamscurrent:
            self._dparamscurrent = value
            self._updateInternals()

    @property
    def dtcurrent(self):
        """Are derivatives with respect to branch lengths current?"""
        return self._dtcurrent

    @dtcurrent.setter
    def dtcurrent(self, value):
        """Set value of `dtcurrent`, update derivatives if needed."""
        assert isinstance(value, bool)
        if value and self.dparamscurrent:
            raise RuntimeError("Can't set both dparamscurrent and dtcurrent True")
        if value != self.dtcurrent:
            self._dtcurrent = value
            self._updateInternals()

    @property
    def t(self):
        """Gets array of branch lengths."""
        return self._t.copy()

    @t.setter
    def t(self, value):
        """Set new branch lengths, update likelihood and derivatives."""
        assert (isinstance(value, scipy.ndarray) and (value.dtype ==
                'float') and (value.shape == self.t.shape))
        if (self._t != value).any():
            self._t = value.copy()
            self._updateInternals()

    @property
    def dloglik(self):
        """`dloglik[param]` is derivative of `loglik` by `param`."""
        assert self.dparamscurrent, "dloglik requires paramscurrent == True"
        return self._dloglik

    @property
    def dloglik_dt(self):
        """`dloglik_dt` is derivative of `loglik` with respect to `t`."""
        assert self.dtcurrent, "dloglik_dt requires dtcurrent == True"
        return self._dloglik_dt

    @property
    def dloglikarray(self):
        """Derivative of `loglik` with respect to `paramsarray`."""
        assert self.dparamscurrent, "dloglikarray requires paramscurrent == True"
        nparams = len(self._index_to_param)
        dloglikarray = scipy.ndarray(shape=(nparams,), dtype='float')
        for (i, param) in self._index_to_param.items():
            if isinstance(param, str):
                dloglikarray[i] = self.dloglik[param]
            elif isinstance(param, tuple):
                dloglikarray[i] = self.dloglik[param[0]][param[1]]
        return dloglikarray

    def updateParams(self, newvalues):
        """Update model parameters and re-compute likelihoods.

        This method is the **only** acceptable way to update model
        parameters. The likelihood is re-computed as needed
        by this method.

        Args:
            `newvalues` (dict)
                A dictionary keyed by param name and with value as new
                value to set. Each parameter name must either be a
                valid model parameter (in `model.freeparams`).
        """
        for (param, value) in newvalues.items():
            if param not in self.model.freeparams:
                raise RuntimeError("Can't handle param: {0}".format(
                        param))
        if newvalues:
            self.model.updateParams(newvalues)
            self._updateInternals()
            self._paramsarray = None

    def _updateInternals(self):
        """Update internal attributes related to likelihood.

        Should be called any time branch lengths or model parameters
        are changed.
        """
        rootnode = self.nnodes - 1
        if self._distributionmodel:
            catweights = self.model.catweights
        else:
            catweights = scipy.ones(1, dtype='float')
        # When there are multiple categories, it is acceptable
        # for some (but not all) of them to have underflow at
        # any given site. Note that we still include a check for
        # Underflow by ensuring that none of the site likelihoods is
        # zero.
        undererrstate = 'ignore' if len(catweights) > 1 else 'raise'
        with scipy.errstate(over='raise', under=undererrstate,
                divide='raise', invalid='raise'):
            self.underflowlogscale.fill(0.0)
            self._computePartialLikelihoods()
            sitelik = scipy.zeros(self.nsites, dtype='float')
            assert (self.L[rootnode] >= 0).all(), str(self.L[rootnode])
            for k in self._catindices:
                sitelik += scipy.sum(self._stationarystate(k) *
                        self.L[rootnode][k], axis=1) * catweights[k]
            assert (sitelik > 0).all(), "Underflow:\n{0}\n{1}".format(
                    sitelik, self.underflowlogscale)
            self.siteloglik = scipy.log(sitelik) + self.underflowlogscale
            self.loglik = scipy.sum(self.siteloglik) + self.model.logprior
            if self.dparamscurrent:
                self._dloglik = {}
                for param in self.model.freeparams:
                    if self._distributionmodel and (param in
                            self.model.distributionparams):
                        name = self.model.distributedparam
                        weighted_dk = (self.model.d_distributionparams[param]
                                * catweights)
                    else:
                        name = param
                        weighted_dk = catweights
                    dsiteloglik = 0
                    for k in self._catindices:
                        dsiteloglik += (scipy.sum(
                                self._dstationarystate(k, name) *
                                self.L[rootnode][k] + self.dL[name][rootnode][k] *
                                self._stationarystate(k), axis=-1) *
                                weighted_dk[k])
                    dsiteloglik /= sitelik
                    self._dloglik[param] = (scipy.sum(dsiteloglik, axis=-1)
                            + self.model.dlogprior(param))
            if self.dtcurrent:
                self._dloglik_dt = 0
                dLnroot_dt = scipy.array([self.dL_dt[n2][rootnode] for
                        n2 in sorted(self.dL_dt.keys())])
                for k in self._catindices:
                    if isinstance(k, int):
                        dLnrootk_dt = dLnroot_dt.swapaxes(0, 1)[k]
                    else:
                        assert k == slice(None)
                        dLnrootk_dt = dLnroot_dt
                    self._dloglik_dt += catweights[k] * scipy.sum(
                            self._stationarystate(k) *
                            dLnrootk_dt, axis=-1)
                self._dloglik_dt /= sitelik
                self._dloglik_dt = scipy.sum(self._dloglik_dt, axis=-1)
                assert self._dloglik_dt.shape == self.t.shape

    def _M(self, k, t, tips=None, gaps=None):
        """Returns matrix exponential `M`."""
        if self._distributionmodel:
            return self.model.M(k, t, tips, gaps)
        else:
            return self.model.M(t, tips, gaps)

    def _dM(self, k, t, param, M, tips=None, gaps=None):
        """Returns derivative of matrix exponential."""
        if self._distributionmodel:
            return self.model.dM(k, t, param, M, tips, gaps)
        else:
            return self.model.dM(t, param, M, tips, gaps)

    def _stationarystate(self, k):
        """Returns the stationarystate ."""
        if self._distributionmodel:
            return self.model.stationarystate(k)
        else:
            return self.model.stationarystate

    def _dstationarystate(self, k, param):
        """Returns the dstationarystate ."""
        if self._distributionmodel:
            return self.model.dstationarystate(k, param)
        else:
            return self.model.dstationarystate(param)

    @property
    def _catindices(self):
        """Returns list of indices of categories."""
        if self._distributionmodel:
            return range(self.model.ncats)
        else:
            return [slice(None)]

    @property
    def _paramlist_PartialLikelihoods(self):
        """List of parameters looped over in `_computePartialLikelihoods`."""
        if self._distributionmodel:
            return [param for param in self.model.freeparams +
                    [self.model.distributedparam] if param not in
                    self.model.distributionparams]
        else:
            return self.model.freeparams

    def _computePartialLikelihoods(self):
        """Update `L`, `dL`, `dL_dt`."""
        for n in range(self.ntips, self.nnodes):
            ni = n - self.ntips # internal node number
            nright = self.rdescend[ni]
            nleft = self.ldescend[ni]
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
            self.L[n] = scipy.ndarray(self._Lshape, dtype='float')
            if self.dparamscurrent:
                for param in self._paramlist_PartialLikelihoods:
                    self.dL[param][n] = scipy.ndarray(self._dLshape[param],
                            dtype='float')
            if self.dtcurrent:
                for n2 in self.dL_dt.keys():
                    self.dL_dt[n2][n] = scipy.zeros(self._Lshape, dtype='float')
            for k in self._catindices:
                if istipr:
                    Mright = MLright = self._M(k, tright,
                            self.tips[nright], self.gaps[nright])
                else:
                    Mright = self._M(k, tright)
                    MLright = broadcastMatrixVectorMultiply(Mright,
                            self.L[nright][k])
                if istipl:
                    Mleft = MLleft = self._M(k, tleft,
                            self.tips[nleft], self.gaps[nleft])
                else:
                    Mleft = self._M(k, tleft)
                    MLleft = broadcastMatrixVectorMultiply(Mleft,
                            self.L[nleft][k])
                self.L[n][k] = MLright * MLleft

                if self.dtcurrent:
                    for (tx, Mx, nx, MLxother, istipx) in [
                            (tright, Mright, nright, MLleft, istipr),
                            (tleft, Mleft, nleft, MLright, istipl)]:
                        if istipx:
                            tipsx = self.tips[nx]
                            gapsx = self.gaps[nx]
                        else:
                            tipsx = gapsx = None
                        dM_dt = self._dM(k, tx, 't', Mx, tipsx, gapsx)
                        if istipx:
                            LdM_dt = dM_dt
                        else:
                            LdM_dt = broadcastMatrixVectorMultiply(
                                    dM_dt, self.L[nx][k])
                        self.dL_dt[nx][n][k] = LdM_dt * MLxother
                        for ndx in self.descendants[nx]:
                            self.dL_dt[ndx][n][k] = broadcastMatrixVectorMultiply(
                                    Mx, self.dL_dt[ndx][nx][k]) * MLxother

                if self.dparamscurrent:
                    for param in self._paramlist_PartialLikelihoods:
                        if istipr:
                            dMright = self._dM(k, tright, param, Mright,
                                    self.tips[nright], self.gaps[nright])
                        else:
                            dMright = self._dM(k, tright, param, Mright)
                        if istipl:
                            dMleft = self._dM(k, tleft, param, Mleft,
                                    self.tips[nleft], self.gaps[nleft])
                        else:
                            dMleft = self._dM(k, tleft, param, Mleft)
                        for j in self._sub_index_param(param):
                            if istipr:
                                dMLright = dMright[j]
                                MdLright = 0
                            else:
                                dMLright = broadcastMatrixVectorMultiply(
                                        dMright[j], self.L[nright][k])
                                MdLright = broadcastMatrixVectorMultiply(
                                        Mright, self.dL[param][nright][k][j])
                            if istipl:
                                dMLleft = dMleft[j]
                                MdLleft = 0
                            else:
                                dMLleft = broadcastMatrixVectorMultiply(
                                        dMleft[j], self.L[nleft][k])
                                MdLleft = broadcastMatrixVectorMultiply(
                                        Mleft, self.dL[param][nleft][k][j])
                            self.dL[param][n][k][j] = ((dMLright + MdLright)
                                    * MLleft + MLright * (dMLleft + MdLleft))

            if ni > 0 and ni % self.underflowfreq == 0:
                # rescale by same amount for each category k
                scale = scipy.amax(scipy.array([scipy.amax(self.L[n][k],
                        axis=1) for k in self._catindices]), axis=0)
                assert scale.shape == (self.nsites,)
                self.underflowlogscale += scipy.log(scale)
                for k in self._catindices:
                    self.L[n][k] /= scale[:, scipy.newaxis]
                    if self.dtcurrent:
                        for n2 in self.dL_dt.keys():
                            self.dL_dt[n2][n][k] /= scale[:, scipy.newaxis]
                    if self.dparamscurrent:
                        for param in self._paramlist_PartialLikelihoods:
                            for j in self._sub_index_param(param):
                                self.dL[param][n][k][j] /= scale[:, scipy.newaxis]

            # free unneeded memory by deleting already used values
            for ntodel in [nright, nleft]:
                if ntodel in self.L:
                    del self.L[ntodel]
                if self.dparamscurrent:
                    for param in self._paramlist_PartialLikelihoods:
                        if ntodel in self.dL[param]:
                            del self.dL[param][ntodel]
                if self.dtcurrent:
                    for n2 in self.dL_dt.keys():
                        if ntodel in self.dL_dt[n2]:
                            del self.dL_dt[n2][ntodel]

    def _sub_index_param(self, param):
        """Returns list of sub-indexes for `param`.

        Used in computing partial likelihoods; loop over these indices."""
        if self._distributionmodel and (param ==
                self.model.distributedparam):
            indices = [()]
        else:
            paramvalue = getattr(self.model, param)
            if isinstance(paramvalue, float):
                indices = [()]
            elif (isinstance(paramvalue, scipy.ndarray) and
                    paramvalue.ndim == 1 and paramvalue.shape[0] > 1):
                indices = [(j,) for j in range(len(paramvalue))]
            else:
                raise RuntimeError("Invalid param: {0}, {1}".format(
                        param, paramvalue))
        return indices



if __name__ == '__main__':
    import doctest
    doctest.testmod()
