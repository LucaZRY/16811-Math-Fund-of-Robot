import numpy as np

"""
PA = LDU
Assumptions:
1. A is a square matrix.
2. A is invertible.


"""

def pa_ldu_decomposition(A):
    A = np.array(A, dtype=float)
    m, n = A.shape

    if m != n:
        raise ValueError("Matrix A must be square.")
    
    A_ = A.copy()

    P = np.eye(m)
    L = np.eye(m)

    for k in range(m):
        # if pivot is zero, swap with a row which below has a non-zero in this column
        if A_[k, k] == 0:
            for i in range(k + 1, m):
                if A_[i, k] != 0:
                    # swap rows in A_
                    A_[[k, i], :] = A_[[i, k], :]
                    # swap rows in P
                    P[[k, i], :] = P[[i, k], :]
                    # swap rows in L (only the part before the k-th column)
                    if k >= 1:
                        L[[k, i], :k] = L[[i, k], :k]
                    break
            else:
                raise ValueError("Matrix A is singular and cannot be decomposed.")
            

def row_change(A_, P_, L_, row1, row2, swap):
    if swap:
        
"""
git check
"""