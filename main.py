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


def interior_point(A, b, c):

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

    #A = np.array([[1,1],[2,1]])
    #b = np.array([12,16])
    #c = -np.array([40, 30])
    m = A.shape[0]
    n = A.shape[1]
    U = 6 * float(max(np.max(np.abs(A)), np.max(np.abs(b)), np.max(np.abs(c))))  # al je to samo max A-ja??
    reduce = 10 ** int(np.log10(U))
    if reduce >= 10:
        return interior_point(A/reduce, b/reduce, c/reduce)

    L = 1 / (((m + 1) * U) ** (3 * (m + 1)))

    M = 4 * n * U / L
    W = ((m + 1) * U) ** (m + 1)
    #W = (m * U) ** m
    #L = 1/(W**2 * 2 * n * ((m + 1) * U)**(3 * (m + 1)))
    #M = 4 * n * U / L

    d = b / W
    rho = d - A.dot(np.ones(A.shape[1]))

    A_prime = np.hstack((A, np.atleast_2d(np.zeros(A.shape[0])).T))

    A_prime = np.vstack((A_prime, np.ones(A_prime.shape[1])))
    # A_prime[-1, -1] = 1
    rho2 = np.hstack((np.atleast_2d(rho), np.atleast_2d(np.array(1))))
    A_prime = np.hstack((A_prime, rho2.T))
    # mu = 1 # TODO?
    mu = 2 * math.sqrt((M ** 2 + np.sum(np.power(c, 2))))  # + 1
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
    delta = 1.0  #(8.0 * math.sqrt(n + 2.0))
    tau = delta / (6 * math.sqrt(n + 2))
    # iz zvezka tau = approx 1/sqrt(n)
    tau = 1 / (6 * math.sqrt(n + 2))  # + 2.0)
    tau = 0.22
    #tau = 0.001
    # print(x_prime, s)
    niter = (math.log(L) + math.log(Q) - (math.log(64 * ((n + 2) ** 2)) + math.log(((m + 1) * U) ** (m + 2))) -
             math.log(mu))/math.log(1 - tau)

    for i in range(int(niter)):
        mu_prime = (1 - tau) * mu
        S = np.diag(s)
        #Sinv = np.linalg.inv(S)
        Sinv = np.diag(1/s)
        X = np.diag(x_prime)

        k = np.linalg.inv(A_prime.dot(Sinv).dot(X).dot(A_prime.T))
        k = k.dot(b_prime - mu_prime * A_prime.dot(Sinv).dot(np.ones(S.shape[0])))

        f = -(A_prime.T.dot(k))
        h = mu_prime * Sinv.dot(np.ones(S.shape[0])) - X.dot(Sinv).dot(f) - x_prime
        x_prime = x_prime + h
        y = y + k
        s = s + f
        mu = mu_prime


    #solution
    Bs = np.nonzero(s < L / (4 * (n + 2)))[0]
    Ns = np.nonzero(x_prime < L / (4 * (n + 2)))[0]
    x_prime[Ns] = 0
    s[Bs] = 0
    x_opt = x_prime[:-2] * W
    s_opt = s[:-2] * W
    y_opt = y[:-2] * W

    Bs = np.nonzero(s_opt < L / (4 * (n + 2)))[0]
    Ns = np.nonzero(x_opt < L / (4 * (n + 2)))[0]

    if len(Bs) < m:
        x_opt, s_opt, optimal_y = interior_point(A[:, Bs], b, c[Bs])
    elif len(Bs) == m:
        x_opt[Bs] = np.linalg.inv(A[:, Bs]).dot(b)

    return x_opt, s_opt, y_opt

if __name__ == '__main__':


    A = np.array([[3, 2, 5, 1],
                  [2, 1, 1, 5],
                  [1, 1, 3, 5],
                  [5, 2, 4, 2]])
    A = A[np.array([0,1,3]), :-1]
    b = np.array([55, 26, 57, 25])
    b = b[:-1]
    c = np.array([20, 10, 15, 32])
    c = c[:-1]
    #A = np.array([[2, 2, 3],
    #              [3, 1, 1],
    #              [1, 1, 1]])
    #A = A[:-1, :-1]
    #b = np.array([140, 150, 133])  # , 104, 125])
    #b = b[:-1]
    #c = np.array([100, 80, 2345])
    #c = c[:-1]
    A = np.concatenate((A, np.eye(A.shape[0])), axis=1)
    c = np.concatenate((c, np.zeros(A.shape[0])))

    x, y, s = interior_point(A, b, c)
    print(x)
    print(c.dot(x))
    print(A.dot(x) <= b)
    res = scipy.optimize.linprog(-c, A_ub=A, b_ub=b)
    print(res.fun)
    print([np.format_float_positional(x, precision=3) for x in res.x[:8]])
    print("-------")
    #A = remove_lin_rows(A)
    # quit()
    A = np.array([[18, 48, 5, 1, 5, 0, 0, 8],
                  [2, 11, 3, 13, 3, 0, 15, 1],
                  [0, 5, 3, 10, 3, 100, 30, 1],
                  [77, 270, 60, 140, 61, 880, 330, 32],
                  [18, 48, 5, 1, 5, 0, 0, 8],
                  [2, 11, 3, 13, 3, 0, 15, 1],
                  [0, 5, 3, 10, 3, 100, 30, 1],
                  [77, 270, 60, 140, 61, 880, 330, 32]])
    c = np.array([10, 22, 15, 45, 40, 20, 87, 21])
    A[:4, :] = -A[:4, :]
    b = np.array([250, 50, 50, 2200,
                  370, 170, 90, 2400])
    b[:4] = -b[:4]
    A = np.concatenate((A, np.eye(A.shape[0])), axis=1)
    c = np.concatenate((c, np.zeros(A.shape[0])))
    res = scipy.optimize.linprog(c, A_ub=A, b_ub=b)
    print(res.fun)
    print([np.format_float_positional(x, precision=3) for x in res.x[:8]])

    ids = A.shape[0]
    A = np.concatenate((A, np.eye(ids)), axis=1)
    c = np.concatenate((c, np.zeros(ids)))

    x, y, s = interior_point(A, b, c)
    np.set_printoptions(linewidth=300)
    print(x.dot(c))
    print([np.format_float_positional(h, precision=3) for h in x[:8]])
    print([np.format_float_positional(h, precision=3) for h in A.dot(x)])
    print([np.format_float_positional(h, precision=3) for h in b])


    # print(A, A.shape    )
    #TODO check intersection of three nearest constraints and put the point there
