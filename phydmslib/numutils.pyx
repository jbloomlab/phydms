"""Numerical routines for ``phydmslib``."""


import numpy
import scipy.linalg
cimport numpy


def broadcastMatrixMultiply(numpy.ndarray a, numpy.ndarray b,
        float alpha=1.0):
    """Broadcast square matrix multiplication using `blas` for speed.

    Args:
        `a` (`numpy.ndarray`, shape `(r, n, n)`)
            First matrix.
        `b` (`numpy.ndarray`, shape `(r, n, n)`)
            Second matrix.
        `alpha` (`float`)
            Multiply each entry in product by this (same meaning
            as for BLAS `dgemm`).

    Returns:
        `ab` (`numpy.ndarray`, shape `(r, n, n)`)
            `ab[r]` is the matrix product of `a[r]`
            with `b[r]`.

    >>> r = 3
    >>> n = 2
    >>> a = numpy.arange(r * n * n, dtype='float').reshape(r, n, n)
    >>> b = numpy.arange(r * n * n, dtype='float').reshape(r, n, n)
    >>> for i in range(r):
    ...   b[i][0][0] = 2
    >>> ab = broadcastMatrixMultiply(a, b)
    >>> expected = scipy.matmul(a, b)
    >>> scipy.allclose(expected, ab)
    True
    """
    assert a.dtype == b.dtype == numpy.double
    assert a.ndim == b.ndim == 3
    cdef int r = a.shape[0]
    assert r == b.shape[0]
    cdef int n = a.shape[1]
    assert n == a.shape[2] == b.shape[2] == b.shape[1]
    if not a.flags['C']:
        a = scipy.ascontiguousarray(a)
    assert b.flags['C'], "b not C contiguous"
    cdef numpy.ndarray ab = numpy.ndarray((r, n, n), dtype=numpy.double)
    cdef int i
    for i in range(r):
        ab[i] = scipy.linalg.blas.dgemm(alpha, a[i], b[i])
        # the overwrite_c option currently appears to be broken, so we do this 
    return ab


def broadcastMatrixVectorMultiply(numpy.ndarray m, numpy.ndarray v,
        float alpha=1.0):
    """Broadcast matrix vector multiplication using `blas` for speed.

    Args:
        `m` (`numpy.ndarray`, shape `(d1, d2, d2)`)
            Array of square matrices to multiply.
        `v` (`numpy.ndarray`, shape `(d1, d2)`)
            Array of vectors to multiply.
        `alpha` (`float`)
            Multiply each entry in product by this (same meaning
            as for BLAS `dgemv`).

    Returns:
        `mv` (`numpy.ndarray`, shape `(d1, d2)`)
            `mv[r]` is the matrix-vector product of `m[r]` with
            `v[r]` for 0 <= `r` <= `d1`.

    >>> m = scipy.array([[[1., 2], [3, 4]], [[5, 6], [7, 8]], [[9, 8], [7, 6]]])
    >>> v = scipy.array([[1., 2], [3, 4], [1, 3]])
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
    assert v.dtype == m.dtype == numpy.double
    assert v.ndim == 2, "v.ndim = {0}".format(v.ndim)
    assert m.ndim == 3
    cdef int r = v.shape[0]
    assert r == m.shape[0]
    cdef int n = v.shape[1]
    assert n == m.shape[1] == m.shape[2]
    if not m.flags['C']:
        m = scipy.ascontiguousarray(m)
    assert v.flags['C']
    cdef numpy.ndarray mv = numpy.ndarray((r, n), dtype=numpy.double)
    cdef int i
    for i in range(r):
        mv[i] = scipy.linalg.blas.dgemv(alpha, m[i], v[i])
    return mv
    #return scipy.sum(m * v[:, None, :], axis=2)


def broadcastGetCols(numpy.ndarray m, numpy.ndarray cols):
    """Get specified columns from `ndarray` of square `ndarrays`.

    This functions uses fast `numpy` broadcasting to get the
    columns from a long array of arrays.

    Args:
        `m` (`numpy.ndarray`, shape `(r, n, n)`
            Array of matrices from which we get the columns
        `cols` (`numpy.ndarray`, type `int`, length `r`)
            All entries should be >= 0 and < `n`.

    Returns:
        `mcols` (`numpy.ndarray`, shape `(r, n)`)
            `mcols[r]` is equal to `mcols[r][cols[r]]`

    >>> n = 2
    >>> r = 3
    >>> m = scipy.arange(r * n * n).reshape(r, n, n)
    >>> cols = scipy.random.random_integers(0, n - 1, r)
    >>> expected = scipy.array([m[i][:, cols[i]] for i in range(r)])
    >>> scipy.allclose(expected, broadcastGetCols(m, cols))
    True
    """
    assert m.ndim == 3
    assert cols.ndim == 1
    assert cols.shape[0] == m.shape[0]
    return m[scipy.arange(m.shape[0]), :, cols]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
