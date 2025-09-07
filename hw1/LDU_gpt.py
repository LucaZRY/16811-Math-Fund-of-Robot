#!/usr/bin/env python3
"""
PA = LDU decomposition with partial pivoting.

- L is unit lower triangular (diag = 1)
- D is diagonal (np.ndarray with only diagonal entries possibly nonzero)
- U is unit upper triangular (diag = 1)
- P is a permutation matrix such that  P A = L D U

Works for singular square matrices (zero pivots allowed).
Does NOT invert D anywhere.
"""

from __future__ import annotations
import numpy as np


def pa_ldu_decomposition(A: np.ndarray, tol: float = 1e-12):
    """
    Compute P, L, D, U such that P @ A = L @ D @ U.

    Parameters
    ----------
    A : (n, n) array_like
        Square matrix (may be singular).
    tol : float
        Tolerance used to detect zero pivots.

    Returns
    -------
    P, L, D, U : np.ndarray
        P (n,n), L (n,n, unit lower), D (n,n diagonal), U (n,n unit upper).

    Notes
    -----
    Algorithm:
        1) Perform Gaussian elimination with *partial pivoting* on A_,
           tracking multipliers in L and row swaps in P.
           This produces an upper-triangular A_ (often called U').
        2) Set D = diag(diag(A_)).
        3) Build U row-by-row by normalizing each row of A_ by the
           corresponding diagonal in D (where nonzero). For zero pivots,
           keep U's row as the standard basis row (diag = 1, rest 0).
    """
    A = np.array(A, dtype=float, copy=True)
    n, m = A.shape
    if n != m:
        raise ValueError("This implementation expects a square matrix.")

    # Working copy (will become an upper-triangular matrix U')
    A_ = A.copy()

    # Initialize outputs
    P = np.eye(n)
    L = np.eye(n)

    # --- Gaussian elimination with partial pivoting (build P, L, and A_ as U') ---
    for k in range(n):
        # Choose pivot row by largest absolute value in column k (from row k down)
        pivot_row = k + np.argmax(np.abs(A_[k:, k]))
        pivot_val = A_[pivot_row, k]

        # If pivot is tiny, treat as zero pivot: skip elimination for this column
        if np.abs(pivot_val) < tol:
            # No elimination possible below; continue to next column
            # (Row k will remain as is; the diagonal A_[k,k] ~ 0)
            continue

        # If we need to swap, do it in A_, P, and the left part of L
        if pivot_row != k:
            A_[[k, pivot_row], :] = A_[[pivot_row, k], :]
            P[[k, pivot_row], :] = P[[pivot_row, k], :]
            # Swap the already-built columns of L (only to the left of k)
            if k > 0:
                L[[k, pivot_row], :k] = L[[pivot_row, k], :k]

        # Eliminate entries below the pivot
        for r in range(k + 1, n):
            if np.abs(A_[k, k]) < tol:
                # Defensive; shouldn't happen because we checked pivot_val
                continue
            mult = A_[r, k] / A_[k, k]
            A_[r, :] -= mult * A_[k, :]
            L[r, k] = mult

    # --- Build D from the diagonal of A_ (U') ---
    D = np.zeros_like(A_)
    diagA = np.diag(A_)
    np.fill_diagonal(D, diagA)

    # --- Build unit-upper U by normalizing rows of A_ by D's diagonal ---
    U = np.eye(n)
    for i in range(n):
        di = D[i, i]
        if np.abs(di) > tol:
            U[i, i:] = A_[i, i:] / di
        else:
            # Zero pivot: the i-th row in A_ should be (numerically) zero from i onward.
            # Keep U[i] = e_i^T (diag 1), but sanity-check numerical zeros:
            if not np.allclose(A_[i, i:], 0.0, atol=1e-9):
                # If this triggers, the matrix is very ill-conditioned for this scheme.
                # You could add complete pivoting or rank-revealing strategies.
                raise ValueError(
                    "Zero pivot encountered but trailing row segment is not ~zero. "
                    "Consider stronger pivoting."
                )

    return P, L, D, U


def check_factorization(P: np.ndarray, L: np.ndarray, D: np.ndarray, U: np.ndarray, A: np.ndarray, tol: float = 1e-9):
    """Return True if P@A ≈ L@D@U within tolerance; otherwise False."""
    left = P @ A
    right = L @ D @ U
    return np.allclose(left, right, atol=tol, rtol=0)


def pretty_print(P, L, D, U, name="A"):
    def fmt(M):
        return np.array2string(M, precision=6, floatmode="maxprec", suppress_small=True)
    print(f"\n=== PA = L D U for {name} ===")
    print("P =\n", fmt(P))
    print("L =\n", fmt(L))
    print("D =\n", fmt(D))
    print("U =\n", fmt(U))


if __name__ == "__main__":
    # ----- Test 1: A3 from your message (singular) -----
    A3 = np.array([
        [1.0, 1.0, 1.0],
        [10.0, 2.0, 9.0],
        [8.0, 0.0, 7.0],
    ])
    P3, L3, D3, U3 = pa_ldu_decomposition(A3)
    pretty_print(P3, L3, D3, U3, name="A3")
    print("Check A3:", check_factorization(P3, L3, D3, U3, A3))


    # ----- Test 2: earlier 3x3 full-rank example -----
    # A1 = np.array([
    #     [10.0, -10.0, 0.0],
    #     [ 0.0,  -4.0, 2.0],
    #     [ 2.0,   0.0, -5.0],
    # ])

    A1 = np.array([[5.0, -5.0,  0.0, 0.0, 0.0],
              [ 5.0, 5.0, 5.0, 0.0, 0.0],
              [ 0.0,  -1.0,  4.0, 1.0, 0.0],
              [ 0.0,  4.0,  -1.0, 2.0, 0.0],
              [ 0.0,  0.0,  2.0, 1.0, 0.0]])
    P1, L1, D1, U1 = pa_ldu_decomposition(A1)
    pretty_print(P1, L1, D1, U1, name="A1")
    print("Check A1:", check_factorization(P1, L1, D1, U1, A1))

    # NOTE: This script handles square matrices (n x n). If you later need
    # a rectangular LDU (m x n, m>=n) like your A2 (5x4), that requires a
    # slightly extended shape-aware version (L is m×m, D is m×n diagonal,
    # U is n×n). Happy to provide that variant if you want.
