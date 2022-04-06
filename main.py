import math
import time

import numpy as np
import scipy.optimize
import scipy.optimize as op
import sympy as sym

def remove_lin_rows(B):
    A = op._remove_redundancy._remove_redundancy_svd(B, np.zeros(B.shape[0]))
    #print(A)
    if A[2] > 0:
        print("unfeasible:", A[1:])
        quit()
    return A[0]


def interior_point():
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
    # A = np.array([[2, 2, 3],
    #               [3, 1, 1],
    #               [1, 1, 1]])
    # A = A[:-1, :-1]
    # #A = np.delete(A, 2, axis=1)
    # b = np.array([140, 150, 133])  # , 104, 125])
    # b = b[:-1]
    # c = -np.array([100, 80, 2345])
    # c = c[:-1]

    A = np.array([[3, 2, 5],
                  [2, 1, 1],
                  #[1, 1, 3],
                  [5, 2, 4]])
    b = np.array([55, 26, 57])
    c = np.array([20, 10, 15])
    A = np.array([[18, 48, 5, 1, 5, 0, 0,  8],
                  [2, 11, 3, 13, 3, 0, 15, 1],
                  [0, 5, 3, 10, 3, 100, 30, 1],
                  [77, 270, 60, 140, 61, 880, 330, 32],
                  [18, 48, 5, 1, 5, 0, 0,  8],
                  [2, 11, 3, 13, 3, 0, 15, 1],
                  [0, 5, 3, 10, 3, 100, 30, 1],
                  [77, 270, 60, 140, 61, 880, 330, 32]])
    c = -np.array([10, 22, 15, 45, 40, 20, 87, 21])
    b = np.array([250, 50, 50, 240, 370, 170, 90, 2200])
    slacks = np.concatenate((np.eye(4), -np.eye(4)), axis=0)
    A = np.concatenate((A, slacks), axis=1)
    #A = np.array([[1,1],[2,1]])
    #b = np.array([12,16])
    #c = -np.array([40, 30])
    A = remove_lin_rows(A)
    #A = np.concatenate((A, np.eye(A.shape[1])), axis=1)
    c = np.concatenate((c, np.zeros(4)))
    m = A.shape[0]
    n = A.shape[1]
    U = float(max(np.max(np.abs(A)), np.max(np.abs(b)), np.max(np.abs(c))))  # al je to samo max A-ja??
    L = 1 / ((m + 1) * U ** (3 * (m + 1)))
    M = 4 * n * U / L
    W = ((m + 1) * U) ** (m + 1)
    d = b / W
    rho = d - A.dot(np.ones(A.shape[1]))

    A_prime = np.hstack((A, np.atleast_2d(np.zeros(A.shape[0])).T))

    A_prime = np.vstack((A_prime, np.ones(A_prime.shape[1])))
    # A_prime[-1, -1] = 1
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
    #tau = 0.01
    # print(x_prime, s)
    while mu >= L * Q / (64 * ((n + 2) ** 2) * ((m + 1) * U) ** (m + 2)):
        if (s == 0).any():
            print("AAAAAAAAAAAA")
        mu_prime = (1 - tau) * mu
        S = np.diag(s)
        #Sinv = np.linalg.inv(S)
        Sinv = np.diag(1/s)
        X = np.diag(x_prime)
        #maybe convert to sym.Matrix???
        k = np.array(A_prime.dot(Sinv).dot(X).dot(A_prime.T).inv()).astype(np.float64)
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
    np.set_printoptions(linewidth=300)
    print((A_prime.dot(x_prime) <= b_prime).all())
    # b_solution = np.zeros(b_prime.shape[0])
    x_solution = x_prime * W

    #invs = check_invs(A, b, c, x_solution, y, s, mu)
    #if not invs.all():
    #    print(invs)
    # b_solution = A_prime[:,Bs].dot(x[Bs])
    # print(b_solution)
    res = scipy.optimize.linprog(c, A_ub=A, b_ub=b)
    print(res)
    print(np.linalg.norm(np.abs(x_solution[:-2] - res.x)))
    print("my solution:   ", x_solution[:-2], c.dot(x_solution[:-2]))
    print("scipy solution:", res.x)
    print(A.dot(x_solution[:-2]), "my solution Ax")
    print(A.dot(res.x), "scipy sol Ax")
    print(b, "b")


if __name__ == '__main__':
    interior_point()
    # print(A, A.shape    )
