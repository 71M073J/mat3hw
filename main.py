import math

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
    # positivity condition
    # assume conditions
    # maximise c = 100x1 + 80x2
    # c -> -100, -80 ker minimiziramo ne maximiziramo
    # 2x1 + 2x2 <= 140
    # 3x1 + 1x2 <= 150
    # A =
    # 2, 2
    # 3, 1
    # b =
    # 140, 150
    # Problem Constraints:
    # Ax ≤ b,????
    # Primal problem

    # min ctx + Mxn+2
    # Ax + pxn+2 = d
    # etx + xn+1 + xn+2 = sum(xn+2) = n + 2
    # x ≥ 0,

    # Dual Problem

    # max bty
    # Aty + s = c
    # s >= 0

    A = np.array([[2, 2],
                  [3, 1]])
    b = np.array([140, 150])
    c = np.array([-100, -80])

    A = remove_lin_rows(A)
    m = A.shape[0]
    n = A.shape[1]
    U = float(max(np.max(np.abs(A)), np.max(np.abs(b)), np.max(np.abs(c))))#al je to samo max A-ja??
    L = 1 / ((m + 1) * U ** (3 * (m + 1)))
    M = 4 * n * U / L
    W = ((m + 1) * U) ** (m + 1)
    d = b / W
    rho = d - A.dot(np.ones(A.shape[1]))

    A_prime = np.vstack((A, np.ones(A.shape[1])))
    A_prime = np.hstack((A_prime, np.atleast_2d(np.zeros(A_prime.shape[0])).T))
    A_prime[-1, -1] = 1
    rho2 = np.hstack((np.atleast_2d(rho), np.atleast_2d(np.array(1))))
    A_prime = np.hstack((A_prime, rho2.T))
    # mu = 1 # TODO?
    mu = math.sqrt(4 * (M ** 2 + np.sum(np.power(c, 2))))  # + 1
    x = np.ones(n + 2)  # starting solution
    # x0[-1] = 1
    s = c + mu
    s = np.hstack((np.atleast_2d(s), np.atleast_2d(np.array([mu, M + mu]))))
    y = np.zeros(m + 1)
    y[-1] = -mu
    b_prime = np.concatenate((d, np.atleast_1d(n + 2)))
    asda = 1
    Q = L/(n + 2)
    while mu >= L * Q / (64 * ((n + 2)**2) * ((m + 1) * U)**(m + 2)):

        xnext = x + h
        ynext = y + k
        snext = s + f




if __name__ == '__main__':
    A = np.array([[2, 3, 5],
                  [1, 1, 1],
                  [2, 2, 2],
                  [2, 1, 1]])
    A = np.random.randint(0, 5, (200, 4))
    A = remove_lin_rows(A)

    A = np.array([[2, 2, -140],  # <= 0 constraints
                  [3, 1, -150]])
    interior_point(cost, A)
    # print(A, A.shape    )
