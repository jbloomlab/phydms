"""Module that wraps an interface to ``Bio++`` for Python.

Written by Jesse Bloom.
"""

import re
import os
import tempfile

import Bio.Phylo
import Bio.Seq
import Bio.Alphabet.IUPAC

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map as cpp_map


cdef extern from "BppExtensions/BppTreeLikelihood.h" namespace "bppextensions":
    cdef cppclass BppTreeLikelihood:
        BppTreeLikelihood(vector[string], vector[string], string, string, bint, cpp_map[int, cpp_map[string, double]], cpp_map[string, double], cpp_map[string, double], bint, bint, bint, bint, char, bint, int, int) except +
        long NSeqs() except +
        long NSites() except +
        void NewickTree(string) except +
        double LogLikelihood() except +
        void OptimizeLikelihood() except +
        cpp_map[string, double] ModelParams() except +
        string OptimizationIgnoredParameters() except +
        cpp_map[string, double] StationaryState(long) except +
        cpp_map[string, double] GetPreferences(long) except +
        void SetPreferences(cpp_map[string, double], long) except +


cdef class PyBppTreeLikelihood:
    """Python wrapped *BppTreeLikelihood* for phylogenetic calculations using ``Bio++``.

    This class wraps a ``C++`` interface to the ``Bio++`` package. It defines
    a single object that can be used for calculating likelihoods, optimizing
    model parameters, optimizing branch lengths, and (for some models) optimizing
    tree topology.

    Objects are instantiated like this:

        *bpptl = BppTreeLikelihood(seqnames, seqs, treefile, model, infertopology, fixedmodelparams, initializemodelparams, oldlikelihoodmethod, fixbrlen, addrateparameter, prefsasparams, recursion, useLog, ngammarates, ncats)*

    where:

        * *seqnames* is a list of strings giving sequence names (must be unique).

        * *seqs* is a list of the same length as *seqnames* giving the sequences
          as upper case ``A``, ``C``, ``G``, ``T``, and ``-`` characters. Must 
          be aligned, and so all must be the same length. They should be
          coding sequences with stop codons trimmed (you may get an undiagnosable
          error if this is not true).

        * *treefile* is an existing file holding a tree in Newick format with tips
          that match those in *seqnames*. Note that this Newick tree can **not**
          have named internal nodes or you will get an error.

        * *model* specifies the substitution model. Allowed values:

            - Variants of the *YNGKP* models of Yang et al, Genetics (2000) 
              (http://www.genetics.org/content/155/1/431.full)
              are specified by a string taking the form *YNGKP_XX_YY* where *XX*
              is method for :math:`\omega` (the dN/dS ratio), and *YY* is method
              to get the equilibrium codon frequencies. For *XX*, the options are:

                - *M0* is a single :math:`\omega` fit by maximum likelihood.

                - *M1* is an :math:`\omega < 1` and an :math:`\omega = 1` with 
                  proportions optimized. This is actually the ``M1a`` model
                  in Wong et al, Genetics (2014) 
                  (http://www.genetics.org/content/168/2/1041)
                  rather than the original ``M1`` in Yang et al, 2000.

                - *M2* is an :math:`\omega < 1`, an :math:`\omega = 1` , and
                  an :math:`\omega > 1` with 
                  proportions optimized. This is actually the ``M2a`` model
                  in Wong et al, Genetics (2014) 
                  (http://www.genetics.org/content/168/2/1041)
                  rather than the original ``M2`` in Yang et al, 2000.

                - *M3* is 3 discrete values of :math:`\omega` with weights and 
                  values optimized by maximum likelihood.
               

                - *M7* is 3 discrete values of :math:`\omega` between 0 and 1 drawn from
                  beta distribution with shape optimized by maximum likelihood.

                - *M8* is like *M7* but with a fourth category with :math:`\omega > 1`.
               
              For *YY*, the valid options are:

                - *fitF3X4* : fit 9 nucleotide frequencies according the F3X4 method.

                - *empF3X4* : frequencies empirically determined from the alignment.

              So for instance, *model* might be *YNGKP_M0_empF3X4*.

            - Experimentally informed substitution models are specified by the
              2-tuple *('ExpCM', aaprefs)* where *aaprefs* is a dictionary
              keyed by integers for every codon site for the sequences in
              *seqs* (1, 2, ... numbering) and the values are dictionaries
              keyed by all 20 amino-acids with values the numerical preference
              for that amino acid at that site.

        * *infertopology* is a Boolean switch specifying if we infer tree topology or
          fix topology to that in *treefile*. Must be *False* if using 
          *YNGKP_M?_F3X4* with *?* > 0.

        * *fixedmodelparams* is a dictionary keyed by names of model parameters to
           fix, with float values being what to fix them to.

        * *initializemodelparams* is a dictionary keyed by the names of model
          parameters to initialize but not to fix, with float values being what
          to initialize them to.

        * *oldlikelihoodmethod* is a Boolean switch that specifies that we use the old
          (not the ``NewLikelihood``) methods of ``Bio++``. Only compatible with
          non-partitioned models.

        * *fixbrlen* is a Boolean switch that specifies that we fix the branch lengths to
          those in *treefile*. Can only be used if *infertopology* is *False*.

        * *addrateparameter* is a Boolean switch specifying that we add a 
          parameter that scales the rates. Only makes sense when *fixbrlen*
          is *True*. 

        * Do we define ExpCM models with preferences as parameters (note that they
          are still not optimized). This is necessary if you want to subsequently
          use *SetPreferences*.

        * *recursion* is the ``Bio++`` likelihood recursion, which can be 'S' (simple)
          or 'D' (double).

        * *useLog* is a Boolean switch specifying if we perform likelihood
          calculations using logarithms.

        * *ngammarates* specifies if we use gamma distributed rates. It
          should be an integer >= 1. If it is one, there is just one constant
          rate. If it is >= 1, then we use discrete gamma-distributed rates
          drawn from this many categories, with the shape parameter
          estimated by maximum likelihood.

        * *ncats* is number of beta-distributed categories for *YNGKP_M7*
          and *YNGKP_M8*
    """

    cdef BppTreeLikelihood *thisptr

    cdef list nts, aminoacids, codons

    cdef dict codon_to_aa

    def __cinit__(self, list seqnames, list seqs, str treefile, model, bint infertopology, fixedmodelparams, initializemodelparams, bint oldlikelihoodmethod, bint fixbrlen, bint addrateparameter, bint prefsasparams, str recursion, bint useLog, int ngammarates, int ncats):
        """Initializes new *PyBppTreeLikelihood* object."""
        # 
        # set up codons, amino acids, nts
        self.nts = list(Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters)
        self.aminoacids = list(Bio.Alphabet.IUPAC.IUPACProtein.letters)
        self.codons = []
        for nt1 in self.nts:
            for nt2 in self.nts:
                for nt3 in self.nts:
                    self.codons.append("%s%s%s" % (nt1, nt2, nt3))
        self.codon_to_aa = dict([(codon, str(Bio.Seq.Seq(codon).translate())) for codon in self.codons])
        #
        # some error checking on calling variables
        if addrateparameter and not fixbrlen:
            raise ValueError("It makes no sense to use addrateparameter without fixbrlen")
        assert len(seqnames) == len(seqs) > 0, "seqnames and seqs must be non-empty and specify the same number of entries"
        assert len(set(seqnames)) == len(seqnames), "seqnames has duplicate entries"
        assert all([isinstance(s, str) for s in seqnames]), "seqnames does not specify entirely strings"
        assert all([(isinstance(s, str) and re.search('^[ATCG-]+$', s) and len(s) % 3 == 0) for s in seqs]), "seqnames does not specify entirely strings with a length a multiple of three and composed of A, T, C, G, and -"
        assert all([len(s) == len(seqs[0]) for s in seqs]), "all sequences in seqs must be the same length"
        assert os.path.isfile(treefile), "treefile of %s does not specify an existing file" % treefile
        assert set([clade.name for clade in Bio.Phylo.read(treefile, 'newick').get_terminals()]) == set(seqnames), "treefile and seqnames do not specify the same set of sequence names"
        assert recursion in ['S', 'D'], "recursion must be 'S' or 'D'"
        assert isinstance(fixedmodelparams, dict) and all([isinstance(key, str) for key in fixedmodelparams.keys()]) and all([isinstance(value, (float, int)) for value in fixedmodelparams.values()]), "fixedmodelparams is not a dictionary keyed by strings with float values"
        assert isinstance(initializemodelparams, dict) and all([isinstance(key, str) for key in initializemodelparams.keys()]) and all([isinstance(value, (float, int)) for value in initializemodelparams.values()]), "initializemodelparams is not a dictionary keyed by strings with float values"
        fixedmodelparams = dict([(key, float(value)) for (key, value) in fixedmodelparams.items()]) # convert any int values to floats
        initializemodelparams = dict([(key, float(value)) for (key, value) in initializemodelparams.items()]) # convert any int values to floats
        #
        # now construct object after processing the model
        yngkp_match = re.compile('^YNGKP_M(?P<modelvariant>\d+)_(emp|fit)F3X4$')
        cdef cpp_map[int, cpp_map[string, double]] preferences  # only needs to be filled with values if using ExpCM
        if isinstance(model, str) and yngkp_match.search(model):
            modelvariant = int(yngkp_match.search(model).group('modelvariant'))
            if modelvariant != 0 and infertopology:
                raise ValueError("Cannot infer topology with %s" % model)
        elif isinstance(model, tuple) and len(model) == 2 and model[0] == 'ExpCM':
            assert isinstance(model[1], dict), "Second entry in model tuple not preferences dict"
            sites = model[1].keys()
            assert len(sites) == len(set(sites)) and min(sites) == 1 and max(sites) == len(seqs[0]) // 3, "Invalid sites in preferences: %s" % str(sites)
            for (r, rprefs) in model[1].items():
                assert isinstance(r, int) and isinstance(rprefs, dict), "preferences must by int keys, dict values"
                assert (sum(rprefs.values()) - 1.0) < 1.0e-4, "preferences do not sum to one"
                for codon in self.codons:
                    aa = self.codon_to_aa[codon]
                    if aa in self.aminoacids:
                        preferences[r][codon] = max(1.0e-6, rprefs[aa]) # don't let any preferences be zero
                    elif aa == '*':
                        preferences[r][codon] = 0.0
                    else:
                        raise ValueError("Invalid aa of %s" % aa)
                    assert isinstance(preferences[r][codon], float), "preference not a number for site %d codon %d" % (r, codon)
            model = 'ExpCM'
        else:
            raise ValueError("Invalid model of %s" % model)
        self.thisptr = new BppTreeLikelihood(seqnames, seqs, treefile, model, infertopology, preferences, fixedmodelparams, initializemodelparams, oldlikelihoodmethod, fixbrlen, addrateparameter, prefsasparams, ord(recursion), useLog, ngammarates, ncats)
        if self.thisptr is NULL:
            raise MemoryError("Failed to allocate pointer to BppTreeLikelihood")

    def __dealloc__(self):
        """Deallocate object"""
        if self.thisptr is not NULL:
            del self.thisptr

    def NSeqs(self):
        """Returns number of sequences."""
        return self.thisptr.NSeqs()

    def NSites(self):
        """Returns number of codon sites."""
        return self.thisptr.NSites()

    def NewickTree(self):
        """Returns Newick representation of current tree as string."""
        try:
            (fd, fname) = tempfile.mkstemp()
            f = os.fdopen(fd, 'w')
            f.close() # close, but we still have the path fname
            self.thisptr.NewickTree(fname)
            with open(fname) as f:
                newicktree = f.read()
        finally:
            try:
                f.close()
            except:
                pass
            if os.path.isfile(fname):
                os.remove(fname)
        return newicktree

    def LogLikelihood(self):
        """Returns current log likelihood (does not perform optimization)."""
        return self.thisptr.LogLikelihood()

    def OptimizeLikelihood(self):
        """Optimizes tree likelihood.

        Optimizes all substitution model parameters, all branch lengths, and
        the topology if *infertopology* was *True* when this object
        was created."""
        self.thisptr.OptimizeLikelihood()

    def ModelParams(self, bint optimizedonly, bint clipnames=True):
        """Returns current values of substitution model parameters.

        If *optimizedonly* is *True*, then only returns parameters being optimized by
        maximum likelihood. If *optimizedonly* is *False*, returns all parameters.
        For instance, there will be non-optimized parameters if the equilibrium codon
        frequencies are being fit empirically.

        If *clipnames* is *True*, clip any prefix from the model parameter name.
        
        The return value is the dictionary *modelparams*, with *modelparams[parametername]*
        being the numerical value of parameter *parametername*.
        """
        modelparams = self.thisptr.ModelParams()
        if optimizedonly:
            ignored = self.thisptr.OptimizationIgnoredParameters()
            if len(ignored):
                for ignore in ignored.split(','):
                    ignorematch = re.compile(ignore.replace('*', '.*'))
                    modelparams = dict([(param, value) for (param, value) in modelparams.items() if not ignorematch.search(param)])
        # clip the prefix (typically "YN98." or "ExpCM." from the model parameter)
        clippedmodelparams = {}
        for modelparam in modelparams:
            if clipnames and ('preferences' not in modelparam) and modelparam != 'Gamma.alpha':
                clippedmodelparam = modelparam[modelparam.index('.') + 1 : ]
            else:
                clippedmodelparam = modelparam
            assert clippedmodelparam not in clippedmodelparams, "Duplicated model param %s after clipping" % clippedmodelparam
            clippedmodelparams[clippedmodelparam] = modelparams[modelparam]
        return clippedmodelparams

    def StationaryState(self, int isite):
        """Returns stationary state of substitution model for site *isite*.

        Numbering is 1 <= *isite* <= *NSites()*

        The return value is the dictionary *stationarystate*, with
        *stationarystate[codon]* giving the stationary state frequency
        of the codon indicated by the upper-case string *codon*.
        """
        stationarystate = self.thisptr.StationaryState(isite)
        assert abs(sum(stationarystate.values()) - 1.0) < 1.0e-5, "The sum of the stationary state is not close to one: %g" % sum(stationarystate.values())
        return stationarystate

    def GetPreferences(self, int isite):
        """Returns the site-specific amino-acid preferences for *isite*.

        Returned as dictionary keyed by amino-acid.
        Only will work (not raise an error) if this
        object was created with an *ExpCM* model.

        Numbering is 1 <= *isite* <= *NSites()*
        """
        return self.thisptr.GetPreferences(isite)


    def SetPreferences(self, cpp_map[string, double] aaprefs, long isite):
        """Sets the site-specific amino-acid preferences to *aaprefs* for *isite*.

        Numbering is 1 <= *isite* <= *NSites()*.

        *aaprefs* is keyed by one-letter amino-acid code, value is preference.

        Only works if this object was created using an *ExpCM* model.
        """
        self.thisptr.SetPreferences(aaprefs, isite)

