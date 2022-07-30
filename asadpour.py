from held_karp import held_karp
import numpy as np


def make_matrices_from_file(M, n):
    dimension = n
    c = []
    for i in range(dimension):
        for j in range(dimension):
            c.append(M[i][j])

    b_eq = np.zeros(2 * dimension, dtype=int)
    A_eq = np.zeros((2 * dimension, dimension * dimension), dtype=int)

    for i in range(dimension):
        b_eq[i] = 1
        for j in range(dimension):
            A_eq[i][i * dimension + j] = 1

    for i in range(dimension):
        b_eq[dimension + i] = 1
        for j in range(dimension):
            A_eq[dimension + i][j * dimension + i] = 1

    b_ub = []
    A_ub = []
    b_ub.append(0)
    A_ub.append([])
    for i in range(dimension * dimension):
        A_ub[-1].append(0)

    return A_eq, b_eq, A_ub, b_ub, c, dimension


def asadpour(M, n):
    A_eq, b_eq, A_ub, b_ub, c, sz = make_matrices_from_file(M, n)
    return held_karp(A_eq, b_eq, A_ub, b_ub, c, sz)
