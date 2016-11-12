"""Utilities for ``phydmslib``."""


import math
import scipy


def BenjaminiHochbergCorrection(pvals, fdr):
    """Benjamini-Hochberg procedure to control false discovery rate.

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



def broadcastMatrixVectorMultiply(m, v):
    """Broadcast matrix vector multiplication.

    This function broadcasts matrix vector multiplication using `scipy`,
    following the approach described here:
    http://stackoverflow.com/questions/26849910/numpy-matrix-multiplication-broadcast

    Args:
        `m` (`numpy.ndarray`, shape `(d1, d2, d2)`)
            Array of square matrices to multiply.
        `v` (`numpy.ndarray`, shape `(d1, d2)`)
            Array of vectors to multiply.

    Returns:
        `mv` (`numpy.ndarray`, shape `(d1, d2)`)
            `mv[r]` is the matrix-vector product of `m[r]` with
            `v[r]` for 0 <= `r` <= `d1`.

    Also can broadcast the case where `v` has an extra dimension preceding
    `(d1, d2)`.

    >>> m = scipy.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]], [[9, 8], [7, 6]]])
    >>> v = scipy.array([[1, 2], [3, 4], [1, 3]])
    >>> mv = scipy.ndarray(v.shape, dtype='int')
    >>> for r in range(v.shape[0]):
    ...   for x in range(v.shape[1]):
    ...      mvrx = 0
    ...      for y in range(v.shape[1]):
    ...        mvrx += m[r][x][y] * v[r][y]
    ...      mv[r][x] = mvrx
    >>> mv2 = broadcastMatrixVectorMultiply(m, v)
    >>> mv.shape == mv2.shape
    True
    >>> scipy.allclose(mv, mv2)
    True
    """
    assert len(m.shape) == 3 and (m.shape[1] == m.shape[2])
    if len(v.shape) == 2:
        assert (v.shape[0] == m.shape[0]) and (v.shape[1] == m.shape[1])
        return scipy.sum(m * v[:, None, :], axis=2)
    elif len(v.shape) == 3:
        assert (v.shape[1] == m.shape[0]) and (v.shape[2] == m.shape[1])
        mv = []
        for i in range(v.shape[0]):
            mv.append(scipy.sum(m * v[i][:, None, :], axis=2))
        return scipy.array(mv)


def broadcastGetCols(m, cols):
    """Get specified columns from `ndarray` of square `ndarrays`.

    This functions uses fast `numpy` broadcasting to get the
    columns from a long array of arrays.

    Args:
        `m` (`numpy.ndarray`, shape `(r, n, n)`
        `cols` (`numpy.ndarray`, type `int`, length `r`)
            All entries should be >= 0 and < `n`.

    Returns:
        `mcols` (`numpy.ndarray`, shape `(r, n)`)
            `mcols[r]` is equal to `mcols[r][cols[r]]`

    Also can broadcast the case where `m` has an extra
    dimension preceding the `(r, n, n)`.

    >>> n = 2
    >>> r = 3
    >>> m = scipy.arange(r * n * n).reshape(r, n, n)
    >>> cols = scipy.random.random_integers(0, n - 1, r)
    >>> expected = scipy.array([m[i][:, cols[i]] for i in range(r)])
    >>> scipy.allclose(expected, broadcastGetCols(m, cols))
    True

    >>> d = 2
    >>> n = 2
    >>> r = 3
    >>> m = scipy.arange(d * r * n * n).reshape(d, r, n, n)
    >>> cols = scipy.random.random_integers(0, n - 1, r)
    >>> expected = []
    >>> for i in range(d):
    ...   expected.append([m[i][j][:, cols[j]] for j in range(r)])
    >>> expected = scipy.array(expected)
    >>> expected.shape == (d, r, n)
    True
    >>> actual = broadcastGetCols(m, cols)
    >>> actual.shape == (d, r, n)
    True
    >>> scipy.allclose(expected, actual)
    True
    """
    assert cols.dtype == 'int'
    if len(m.shape) == 3:
        (r, nx, ny) = m.shape
        assert nx == ny
        assert cols.shape == (r,)
        return m[scipy.arange(r), :, cols]
    else:
        (d, r, nx, ny) = m.shape
        assert nx == ny
        assert cols.shape == (r,)
        mcols = []
        for i in range(d):
            mcols.append(m[i][scipy.arange(r), :, cols])
        return scipy.array(mcols)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
