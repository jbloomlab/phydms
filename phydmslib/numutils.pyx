"""Numerical routines for ``phydmslib``."""


import scipy


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
    else:
        raise RuntimeError("invalid shape")


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
