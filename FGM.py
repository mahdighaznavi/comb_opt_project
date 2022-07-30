# %%
import numpy as np
from scipy import optimize
from constants import INF


# %%
def make_matrices_from_file(file_path):
    a = open(file_path)
    a.readline()
    a.readline()
    a.readline()
    dimension = int(a.readline().split()[1])
    while a.readline().strip() != "EDGE_WEIGHT_SECTION":
        continue
    e = []
    for p in range(dimension):
        e.append([])
        while len(e[-1]) < dimension:
            tmp = map(int, a.readline().split())
            for k in tmp:
                if k == "100000000":
                    k = "9999999"
                e[-1].append(k)

    return dimension, e


# %%
# MMAX = 10000
# INF = 9999999
# index_Edges = {}
# n, M = make_matrices_from_file('ftv33.atsp')
# m = n ** 2
#

# %%
def minimum_Matching(graph, sz):
    row_ind, col_ind = optimize.linear_sum_assignment(graph)
    nodes = []
    mark = np.zeros(sz, dtype=int)
    comp = np.zeros(sz, dtype=int)
    x = 0
    for i in range(sz):
        t = i
        if mark[t] == 0:
            x += 1
        while not mark[t]:
            comp[t] = x - 1
            nodes.append(t)
            mark[t] = 1
            t = col_ind[t] - sz
    cost = graph[row_ind, col_ind].sum()
    return cost, nodes, comp


# %%
def split_Graph(graph, sz):
    new_graph = np.zeros((2 * sz, 2 * sz), dtype=int)
    for i in range(sz):
        for j in range(sz):
            new_graph[i][j] = INF
    for i in range(sz):
        for j in range(sz):
            new_graph[i][sz + j] = graph[i][j]
    return new_graph


# %%
def contract(graph, comp, sz):
    new_size = 0
    for i in range(sz):
        new_size = max(new_size, comp[i])
    new_size += 1
    new_g = np.zeros((new_size, new_size), dtype=int)
    for i in range(new_size):
        new_g[i][i] = INF
    in_edges = np.zeros(new_size, dtype=int)
    out_edges = np.zeros(new_size, dtype=int)
    for i in range(sz):
        for j in range(sz):
            if comp[i] != comp[j]:
                new_g[comp[i]][comp[j]] = max(new_g[comp[i]][comp[j]], graph[i][j])
    for i in range(sz):
        for j in range(sz):
            if comp[j] - comp[i] == 1:
                if new_g[comp[i]][comp[j]] == graph[i][j]:
                    out_edges[comp[i]] = i
                    in_edges[comp[j]] = j
            if (comp[i] == new_size - 1) and (comp[j] == 0):
                out_edges[comp[i]] = i
                in_edges[comp[j]] = j
    added = 0
    for i in range(new_size):
        if in_edges[i] == out_edges[i]:
            continue
        added += graph[in_edges[i]][out_edges[i]]
    return new_g, new_size, added


# %%
def solve(graph, sz):
    cost, nodes, comp = minimum_Matching(split_Graph(graph, sz), sz)
    new_graph, new_size, added = contract(graph, comp, sz)
    if new_size == 1:
        return cost
    return cost + added + solve(new_graph, new_size)

