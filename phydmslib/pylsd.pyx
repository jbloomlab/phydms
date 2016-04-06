"""Cython module that wraps a Python interface to ``LSD``.

``LSD`` is a least-squares dating method for phylogenies.
See https://dx.doi.org/10.1093/sysbio/syv068 for
the reference.

This module wraps the ``LSD`` program distributed with ``phydms``.

Written by Jesse Bloom.
"""

import re
import os
import tempfile
import Bio.Phylo


cdef extern from "LSDExtensions/lsd.h":
    int lsd(int, char* []) except +


def RunLSD(intree, dates, seqlength, outtree):
    """Dates sequences in phylogeny using least-square criteria.

    This function does the equivalent of running the ``lsd`` 
    program as wrapped by this module. It gives essentially the
    same result you would expect from running::

        lsd -i intree -d dates -v -s seqlength -c -r as -o outtree

    Essentially, if each sequence has a date assigned to it, this
    function finds the root and dates internal nodes.

    The calling arguments are as follows:

    * *intree* is the name of an existing file containing the 
      input tree in Newick format.
      Spaces are not allowed in the sequence names.

    * *dates* is a dictionary giving the dates (as numbers) for each tip
      in *intree*. For each of these tips, there should be a key in
      *dates* that has as its value the numeric date for that sequence.
      There shouldn't be any extraneous keys.

    * *seqlength* is the length of the sequence.

    * *outtree* is a string giving the name of the file that we create
      that contains the dated tree. You will get an error if this file
      already exists, so remove any existing files of this name.
      In this file, the branch lengths are now in units of elapsed time
      that match those provided in *datefile*. The is the ``.date.newick``
      file produced by ``lsd``.

    The return value of this function is a number giving the date of the
    most recent common ancestor.
    """
    assert not os.path.isfile(outtree), "outtree %s already exists. You must delete this file before calling RunLSD." % outtree
    assert isinstance(seqlength, int) and 1 <= seqlength, "seqlength must be integer >= 1"
    assert os.path.isfile(intree), "intree %s does not exist" % intree

    tipnames = [clade.name for clade in Bio.Phylo.read(intree, 'newick').get_terminals()]
    assert len(set(tipnames)) == len(tipnames), "the tips in %s have duplicate names" % intree
    tipnames = set(tipnames)
    abberrantnames = map(str, tipnames.symmetric_difference(set(dates.keys())))
    if abberrantnames:
        raise ValueError("intree %s and dates do not specify exactly matching sets of sequence names.\nThe mismatched ones are:\n%s" % (intree, ', '.join(abberrantnames)))
    assert all([not re.search('\s', tip) for tip in tipnames]), "Some of the tips specified in intree and dates have spaces -- this is not allowed"

    outtreedir = tempfile.mkdtemp()
    outtreeresults = "%s/outtree.results"  % outtreedir # name of file containing results
    outtreedated = '%s/outtree.results.date.newick'  % outtreedir # name of file containing dated tree
    (fd_dates, tempdatesfile) = tempfile.mkstemp()
    with os.fdopen(fd_dates, 'w') as f:
        f.write('%d\n%s' % (len(dates), '\n'.join(['%s\t%g' % tup for tup in dates.items()])))
    cdef char* c_argv[13] # bit of a hack to hardcode number of arguments
    c_argv[0] = b'lsd'
    c_argv[1] = b'-i'
    arg2 = bytes(intree)
    c_argv[2] = arg2
    c_argv[3] = b'-d'
    arg4 = bytes(tempdatesfile)
    c_argv[4] = arg4
    c_argv[5] = b'-v'
    c_argv[6] = b'-s'
    arg7 = bytes(str(seqlength))
    c_argv[7] = arg7
    c_argv[8] = b'-c'
    c_argv[9] = b'-r'
    c_argv[10] = b'as'
    c_argv[11] = b'-o'
    arg12 = bytes(outtreeresults)
    c_argv[12] = arg12
    try:
        returncode = lsd(13, c_argv) # now using hardcoded number of arguments
        if returncode != 0:
            raise RuntimeError("Running ``lsd`` via the cython wrapper gave a non-zero return code.")
        assert os.path.isfile(outtreeresults) and os.path.isfile(outtreedated), "Failed to find expected files %s and %s\nFiles in outtreedir %s are:\n%s" % (outtreeresults, outtreedated, outtreedir, ', '.join(os.listdir(outtreedir)))
        with open(outtreedated) as f:
            treetext = f.read()
        tmrca_match = re.compile('Tree 1 rate\s(?P<rate>\d+\.\d+),{0,1}\stMRCA\s(?P<tmrca>\-{0,1}\d+\.\d+),{0,1}\sobjective[_\s]function\s(?P<objective>\d+\.\d+)')
        with open(outtreeresults) as f:
            resultstext = f.read()
        m = tmrca_match.search(resultstext)
        assert m, "Failed to match tMRCA in lsd results. Results are:\n%s" % resultstext
        tmrca = float(m.group('tmrca'))
    finally:
        if os.path.isfile(tempdatesfile):
            os.remove(tempdatesfile)
        if os.path.isdir(outtreedir):
            for fname in os.listdir(outtreedir):
                os.remove("%s/%s" % (outtreedir, fname))
            os.rmdir(outtreedir)

    with open(outtree, 'w') as f:
        f.write(treetext)
    assert tipnames == set([clade.name for clade in Bio.Phylo.read(outtree, 'newick').get_terminals()]), "outtree %s does not have all of the expected tip nodes" % outtree

    return tmrca


if __name__ == '__main__':
    import doctest
    doctest.testmod()
