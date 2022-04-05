import math
import time

import numpy as np
import scipy.optimize
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
    # solution: x1 = 40, x2 = 30
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
    # introduce another unknown variable into A?
    A = np.array([[2, 2, 3, 1, 0],
                  [3, 1, 1, 0, 1],
                  [1, 1, 1, 0, 0]])
    b = np.array([140, 150, 133])#, 104, 125])
    c = -np.array([100, 80, 2345, 0, 0])

    A = remove_lin_rows(A)
    m = A.shape[0]
    n = A.shape[1]
    U = float(max(np.max(np.abs(A)), np.max(np.abs(b)), np.max(np.abs(c))))  # al je to samo max A-ja??
    L = 1 / ((m + 1) * U ** (3 * (m + 1)))
    M = 4 * n * U / L
    W = ((m + 1) * U) ** (m + 1)
    d = b / W
    print(d)
    print(np.ones(A.shape[0]).dot(A))
    rho = d - A.dot(np.ones(A.shape[1]))

    A_prime = np.vstack((A, np.ones(A.shape[1])))
    A_prime = np.hstack((A_prime, np.atleast_2d(np.zeros(A_prime.shape[0])).T))
    A_prime[-1, -1] = 1
    rho2 = np.hstack((np.atleast_2d(rho), np.atleast_2d(np.array(1))))
    A_prime = np.hstack((A_prime, rho2.T))
    # mu = 1 # TODO?
    mu = math.sqrt(4 * (M ** 2 + np.sum(np.power(c, 2))))  # + 1
    x_prime = np.ones(n + 2)  # starting solution
    # x0[-1] = 1
    s = c + mu
    s = np.hstack((np.atleast_2d(s), np.atleast_2d(np.array([mu, M + mu])))).flatten()
    y = np.zeros(m + 1)
    y[-1] = -mu
    b_prime = np.concatenate((d, np.atleast_1d(n + 2)))
    Q = L / (n + 2)
    tempmu = mu
    # za delta <= 1/6
    delta = 1.0 / (8.0 * math.sqrt(n + 2.0))
    # tau = delta / math.sqrt(n)
    # iz zvezka tau = approx 1/sqrt(n)
    tau = 1 / math.sqrt(n + 2)  # + 2.0)
    # print(x_prime, s)
    while mu >= L * Q / (64 * ((n + 2) ** 2) * ((m + 1) * U) ** (m + 2)):
        mu_prime = (1 - tau) * mu
        S = np.diag(s)
        Sinv = np.linalg.inv(S)
        X = np.diag(x_prime)
        k = np.linalg.inv(A_prime.dot(Sinv).dot(X).dot(A_prime.T))
        k = k.dot(b_prime - mu_prime * A_prime.dot(Sinv).dot(np.ones(S.shape[0])))
        f = -(A_prime.T.dot(k))
        h = mu_prime * Sinv.dot(np.ones(S.shape[0])) - X.dot(Sinv).dot(f) - x_prime
        x_prime = x_prime + h
        y = y + k
        s = s + f
        mu = mu_prime

    # print(A_prime.dot(x), b_prime)
    # Bs = np.nonzero(s < L / (4 * (n + 2)))[0]
    # Ns = np.nonzero(x_prime < L / (4 * (n + 2)))[0]
    # print(Bs, Ns)
    print((A_prime.dot(x_prime) <= b_prime).all())
    # b_solution = np.zeros(b_prime.shape[0])
    x_solution = x_prime * W
    print(x_solution)
    # b_solution = A_prime[:,Bs].dot(x[Bs])
    # print(b_solution)
    for i in range(A.shape[0]):
        sumcost = 0
        for k in range(A.shape[1]):
            sumcost += A[i, k] * x_solution[k]
        print(sumcost, b[i])
    print(x_solution[:-2])
    print(c.dot(x_solution[:-2]))

    #res = scipy.optimize.linprog(c, A_ub=A[:-1], b_ub=b[:-1], A_eq=np.atleast_2d(A[-1]), b_eq=np.atleast_1d(b[-1]))
    #print(c.dot(res.x))
    #print(np.linalg.norm(np.abs(x_solution[:-2] - res.x)) < 1e-5)



if __name__ == '__main__':
    interior_point()
    # print(A, A.shape    )
