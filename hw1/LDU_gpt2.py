#!/usr/bin/env python3
"""
Rectangular PA = L D U decomposition (m x n, m >= n).

- P : (m x m) permutation matrix
- L : (m x m) unit lower triangular
- D : (m x n) "diagonal" rectangular matrix whose nonzero entries are only on
      the main diagonal min(m,n); i.e., D[i,i] are the pivots, others zero
- U : (n x n) unit upper triangular

Satisfies:  P @ A = L @ D @ U

Notes:
- Works with singular problems (zero pivots allowed).
- Uses partial pivoting (by column).
"""

from __future__ import annotations
import numpy as np


def pa_ldu_rectangular(A: np.ndarray, tol: float = 1e-12):
    """
    Compute P, L, D, U for a rectangular A (m x n) with m >= n such that:
        P @ A = L @ D @ U

    Parameters
    ----------
    A : array_like, shape (m, n)
    tol : float
        Threshold to treat a pivot as zero.

    Returns
    -------
    P : (m, m) ndarray
    L : (m, m) ndarray (unit lower)
    D : (m, n) ndarray (rectangular diagonal)
    U : (n, n) ndarray (unit upper)
    """
    A = np.array(A, dtype=float, copy=True)
    m, n = A.shape
    if m < n:
        raise ValueError("Require m >= n for this rectangular LDU.")

    # Working copy; we will do elimination on rows (like producing an upper-trapezoidal R)
    R = A.copy()

    P = np.eye(m)
    L = np.eye(m)

    # Gaussian elimination with partial pivoting for each column k=0..n-1
    for k in range(n):
        # Choose the pivot row among rows k..m-1 in column k
        pivot_row = k + np.argmax(np.abs(R[k:, k]))
        pivot_val = R[pivot_row, k]

        # If pivot is tiny, treat as zero -> skip elimination for this column
        if abs(pivot_val) < tol:
            # leave R[k:,k] as-is; D[k,k] will be zero later
            continue

        # Swap pivot row into position k in R, and mirror swaps in P and L (to the left of k)
        if pivot_row != k:
            R[[k, pivot_row], :] = R[[pivot_row, k], :]
            P[[k, pivot_row], :] = P[[pivot_row, k], :]
            if k > 0:
                L[[k, pivot_row], :k] = L[[pivot_row, k], :k]

        # Eliminate entries below the pivot (rows k+1..m-1)
        for r in range(k + 1, m):
            if abs(R[k, k]) < tol:
                continue
            mult = R[r, k] / R[k, k]
            R[r, :] -= mult * R[k, :]
            L[r, k] = mult

    # Build rectangular diagonal D (m x n) with D[i,i] = R[i,i] for i<min(m,n)=n
    D = np.zeros((m, n), dtype=float)
    for i in range(n):
        D[i, i] = R[i, i]

    # Build unit-upper U (n x n) by normalizing the first n rows of R by D[i,i]
    U = np.eye(n)
    for i in range(n):
        di = D[i, i]
        if abs(di) > tol:
            U[i, i:] = R[i, i:] / di
        else:
            # Zero pivot -> keep U[i] = e_i^T; sanity check trailing segment ~ 0
            if not np.allclose(R[i, i:], 0.0, atol=1e-9):
                raise ValueError(
                    "Zero pivot with nonzero trailing row; consider stronger pivoting."
                )

    return P, L, D, U


def check_rect(P, L, D, U, A, tol: float = 1e-9):
    """Check P@A ≈ L@D@U for rectangular case."""
    return np.allclose(P @ A, L @ D @ U, atol=tol, rtol=0)


def pretty(M, name):
    print(f"{name} =\n{np.array2string(M, precision=6, floatmode='maxprec', suppress_small=True)}\n")


if __name__ == "__main__":
    # ---- Example: your A2 (5x4) ----
    A2 = np.array([
        [5., -5.,  0., 0.],
        [5.,  5.,  5., 0.],
        [0., -1.,  4., 1.],
        [0.,  4., -1., 2.],
        [0.,  0.,  2., 1.],
    ])

    P2, L2, D2, U2 = pa_ldu_rectangular(A2)
    pretty(P2, "P2")
    pretty(L2, "L2")
    pretty(D2, "D2")
    pretty(U2, "U2")
    print("Check A2:", check_rect(P2, L2, D2, U2, A2))

    # You can also test it on square matrices (m=n); it coincides with the square algorithm.
    A1 = np.array([
        [10., -10.,  0.],
        [ 0.,  -4.,  2.],
        [ 2.,   0., -5.],
    ])
    P1, L1, D1, U1 = pa_ldu_rectangular(A1)
    pretty(P1, "P1")
    pretty(L1, "L1")
    pretty(D1, "D1")
    pretty(U1, "U1")
    print("Check A1:", check_rect(P1, L1, D1, U1, A1))
