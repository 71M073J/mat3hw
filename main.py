import numpy as np
import scipy.optimize as op


def remove_lin_rows(B):
    A = op._remove_redundancy._remove_redundancy_svd(B, np.zeros(B.shape[0]))
    if A[2] > 0: print("infeasible:", A[1:])
    return A[0]


#     A = B.astype(float)
#     for allRows in range(A.shape[1]):  # kateri row odštevamo
#         for kateriRow in range(allRows + 1, A.shape[0]):
#             if A[allRows, allRows] != 0:
#                 A[kateriRow, :] -= (A[kateriRow, allRows] / A[allRows, allRows]) * A[allRows, :]
#     indexes = A.any(axis=1)
#     return B[indexes]


def cost():
    ...


def interior_point(cost=cost, A=None):
    # x is the solution vector for which we are solving
    # c is the cost function
    # Ax = b
    # A and b are the equality conditions
    # x >= 0
    #positivity condition
    #assume conditions
    # c = 100x1 + 80x2
    #2x1 + 2x2 <= 140
    #3x1 + 1x2 <= 150

    #Problem Constraints:
    #Ax ≤ b,????

    #Primal problem

    #min ctx + Mxn+2
    #Ax + pxn+2 = d
    #etx + xn+1 + xn+2 = sum(xn+2) = n + 2
    #x ≥ 0,

    #Dual Problem

    # max bty
    #Aty + s = c
    #s >= 0




    A = remove_lin_rows(A)
    m = A.shape[0]
    n = A.shape[1]
    U = max(A)

    M = 4 * n * U / L
    x0 = np.zeros(A.shape[1] + 2)  # starting solution
    x0[-1] = 1
    c = np.zeros(A.shape[1])


if __name__ == '__main__':
    A = np.array([[2, 3, 5],
                  [1, 1, 1],
                  [2, 2, 2],
                  [2, 1, 1]])
    A = np.random.randint(0, 5, (200, 4))
    A = remove_lin_rows(A)

    A = np.array([[2,2,-140], # <= 0 constraints
                  [3,1,-150]])
    interior_point(cost, A)
    # print(A, A.shape    )
