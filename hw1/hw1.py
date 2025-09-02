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
        if A_[k, k] == 0:
            row_to_change = k
            for t in range(k + 1, m):
                if A_[t, k] != 0:
                    row_to_change = t
                    break
            P, A_, L = row_change(A_, P, L, row_to_change, k, is_swap=True)

        for t in range(k + 1, m):
            P, A_, L = row_change(A_, P, L, t, k, is_swap=False)

    # Extract D and U such that A_ = D @ U
    D = np.diag(np.diag(A_))
    U = np.linalg.inv(D) @ A_

    return P, L, D, U

def row_change(A_, P_, L_, row2, row1, is_swap):
    """
    Handle row swapping or Gaussian elimination for LDU.
    Args:
        A_, P_, L_ (ndarray): Current matrices
        row2, row1 (int): Row indices
        is_swap (bool): Whether to perform row swap or elimination
    Returns:
        Updated (P_, A_, L_)
    """
    A_ = A_.copy()
    P_ = P_.copy()
    L_ = L_.copy()

    if is_swap:
        A_[[row1, row2], :] = A_[[row2, row1], :]
        P_[[row1, row2], :] = P_[[row2, row1], :]
        L_[[row1, row2], :row1] = L_[[row2, row1], :row1]
    else:
        times = A_[row2, row1] / A_[row1, row1]
        A_[row2, :] -= times * A_[row1, :]
        L_[row2, row1] = times

    return P_, A_, L_


# Test Case:
# A = np.array([[2.0, 1.0, 1.0],
#               [4.0, -6.0, 0.0],
#               [-2.0, 7.0, 2.0]])

# A = np.array([[ 7.0, -1.0,  4.0],
#               [ 2.0, -8.0, -3.0],
#               [-2.0,  6.0, -5.0]])


# A = np.array([[1.0, 1.0, 1.0],
#               [0.0, 2.0, 5.0],
#               [2.0, 5.0, -1.0]])

A = np.array([[3.0, 1.0, 2.0],
              [0.0, 1.0, 0.0],
              [0.0, 0.0, 1.0]])

# A = np.array([[-9.0, -3.0,  4.0],
#               [ 0.0, -6.0, -2.0],
#               [ 2.0,  3.0,  5.0]])



P, L, D, U = pa_ldu_decomposition(A)

print("P =\n", P)
print("L =\n", L)
print("D =\n", D)
print("U =\n", U)

PA = P @ A
LDU = L @ D @ U
print("PA ≈ LDU:", np.allclose(PA, LDU))
