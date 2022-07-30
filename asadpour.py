from held_karp import held_karp


def make_matrices_from_file(file_path):
    a = open(file_path)
    a.readline()
    a.readline()
    a.readline()
    dimension = int(a.readline().split()[1])
    while a.readline().strip() != "EDGE_WEIGHT_SECTION":
        continue
    e = []
    c = []
    for i in range(dimension):
        e.append([])
        while len(e[-1]) < dimension:
            tmp = map(int, a.readline().split())
            for t in tmp:
                c.append(t)
                e[-1].append(t)

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

    return A_eq, b_eq, A_ub, b_ub, c, dimension, e


def asadpour(file_path):
    A_eq, b_eq, A_ub, b_ub, c, sz, M = make_matrices_from_file(file_path)
    held_karp(A_eq, b_eq, A_ub, b_ub, c, sz)
