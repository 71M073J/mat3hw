import numpy as np


def remove_lin_rows(A):
    return A


def interior_point(A):
    A = remove_lin_rows(A)
    m = A.shape[0]
    n = A.shape[1]
    U = max(A)
    
    M = 4 * n * U / L
    x0 = np.zeros(A.shape[1] + 2) #starting solution
    x0[-1] = 1
    c = np.zeros(A.shape[1])


if __name__ == '__main__':
    A = np.array([[2, 3, 5],
                  [1, 1, 1],
                  [2, 1, 1]])
