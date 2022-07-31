import maxflow
import numpy as np
from scipy import optimize
from constants import MMAX, INF


# %%
def solve_LP(c, A_ub, b_ub, A_eq, b_eq):
    return optimize.linprog(c=c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, method='interior-point', )


# %%
def add_violated_constraint(nodes, size):
    constraint = []
    for i in range(size ** 2):
        u = i // size
        v = i % size
        if (u in nodes) and (not (v in nodes)):
            constraint.append(-1)
        else:
            constraint.append(0)
    return constraint


# %%
def contract(graph, nodes, size):
    size = size - len(nodes) + 1
    ind = np.zeros(size, dtype=int)
    ctr_graph = np.zeros((size, size), dtype=int)
    t = 0
    for i in range(size):
        if i in nodes:
            ind[i] = size - 1
        else:
            ind[i] = t
            t += 1
    for i in range(size):
        for j in range(size):
            if (i in nodes) and (not (j in nodes)):
                ctr_graph[size - 1][ind[j]] += graph[i][j]
            if (not (i in nodes)) and (j in nodes):
                ctr_graph[ind[i]][size - 1] += graph[i][j]
            if (not (i in nodes)) and (not (j in nodes)):
                ctr_graph[ind[i]][ind[j]] += graph[i][j]
    return ctr_graph


def create_maxflow_graph(graph, s, t):
    n = len(graph)
    g = maxflow.Graph[float](n - 2, (n - 2) * (n - 2))
    nodes = g.add_nodes(n - 2)
    real_i = 0
    for i in range(n):
        if i == s or i == t:
            continue
        real_j = 0
        for j in range(n):
            if j == s or j == t:
                continue
            g.add_edge(nodes[real_i], nodes[real_j], graph[i][j], graph[j][i])
            real_j += 1

        real_i += 1

    real_i = 0
    for i in range(n):
        if i == s or i == t:
            continue
        g.add_tedge(nodes[real_i], graph[s][i], graph[i][t])
        real_i += 1
    return g, nodes


# %%
def get_min_cut(graph, s, t, return_nodes=False):
    g, nodes = create_maxflow_graph(graph, s, t)
    flow = g.maxflow() + graph[s][t]
    if return_nodes:
        cut = [s]
        real_i = 0
        for i in range(len(graph)):
            if i == s or i == t:
                continue
            if g.get_segment(nodes[real_i]) == 0:
                cut.append(i)
            real_i += 1
        return flow, cut
    return flow


# %%
def min_cut(graph, size):
    s = 0
    temp = INF
    for j in range(size):
        if j == s:
            continue
        temp = min(temp, get_min_cut(graph, s, j))
        temp = min(temp, get_min_cut(graph, j, s))
    return temp


# %%
def find_violated_cut(graph, size):
    s = t = 0
    nodes = []
    f = min_cut(graph, size)
    if f < MMAX:
        for i in range(size):
            if i == s:
                continue
            if get_min_cut(graph, s, i) <= f:
                t = i
            if get_min_cut(graph, i, s) <= f:
                t = s
                s = i
        flow, nodes = get_min_cut(graph, s, t, True)
    return nodes


# %%
def make_graph(edges, size):
    g = []
    for i in range(size):
        g.append([])
        for j in range(size):
            g[-1].append(0)
    for i in range(size ** 2):
        g[i // size][i % size] = edges[i]
    return g


def held_karp(A_eq, b_eq, A_ub, b_ub, c, sz, return_complete=False):
    x_star = solve_LP(c=c, A_eq=A_eq, b_eq=b_eq, A_ub=A_ub, b_ub=b_ub)
    mg = make_graph(x_star.x, sz)
    while min_cut(mg, sz) < MMAX:
        T = find_violated_cut(mg, sz)
        if len(T) < 1:
            break
        b_ub.append(-1)
        A_ub.append(add_violated_constraint(T, sz))
        x_star = solve_LP(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq)
        mg = make_graph(x_star.x, sz)

        print(x_star.fun)
    if return_complete:
        return x_star, A_ub, b_ub
    else:
        return x_star.fun, x_star.x
