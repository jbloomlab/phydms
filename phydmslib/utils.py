"""Utilities for ``phydmslib``."""


import math


def BenjaminiHochbergCorrection(pvals, fdr):
    """Implements Benjamini-Hochberg procedure to control false discovery rate.

    Calling arguments:
 
    *pvals* : a list of tuples of *(label, p)*  where *label* is some label assigned
    to each data point, and *p* is the corresponding *P-value*.

    *fdr* : the desired false discovery rate

    The return value is the 2-tuple *(pcutoff, significantlabels)*. After applying
    the algorithm, all data points with *p <= pcutoff* are declared significant.
    The labels for these data points are in *significantlabels*. If there are no
    significant sites, *pcutoff* is returned as the maximum P-value that would
    have made a single point significant.
    """
    num_tests = len(pvals)

    # sort by p-value
    sorted_tests = sorted(pvals, key=lambda tup: tup[1])

    # find maximum rank for which p <= (rank/num_tests)*FDR
    max_rank = 0
    pcutoff = None
    for (rank, (label, p)) in enumerate(sorted_tests):
        rank = rank + 1 # rank beginning with 1 for smallest p-value (there is no rank 0)
        bh_threshold = fdr * float(rank) / num_tests
        if p <= bh_threshold: 
            assert rank > max_rank
            max_rank = rank
            pcutoff = bh_threshold

    # pcutoff to have one significant site if there are none
    if pcutoff == None:
        pcutoff = 1.0 / num_tests * fdr

    # collect significant ranks:
    significantlabels = []
    for (rank, (label, p)) in enumerate(sorted_tests):
        rank = rank + 1 # rank beginning with 1 for site with smallest p-vaalue
        if rank <= max_rank:
            assert p <= pcutoff
            significantlabels.append(label)

    return (pcutoff, significantlabels)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
