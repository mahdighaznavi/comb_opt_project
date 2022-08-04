import math
import random
from math import exp
from math import log as ln

import networkx as nx
import numpy as np

from constants import INF
from held_karp import held_karp


def dfs(graph, mark, x, v, sz):
    mark[v] = x
    for i in range(sz):
        if i == v:
            continue
        if (not mark[i]) and graph[v][i] == 1:
            dfs(graph, mark, x, i, sz)


# %%
def make_components(graph, comp, sz):
    new_sz = 0
    for i in range(sz):
        new_sz = max(new_sz, comp[i])
    new_sz += 1
    new_g = np.zeros((new_sz, new_sz), dtype=int)
    for i in range(sz):
        for j in range(sz):
            if comp[i] != comp[j]:
                new_g[comp[i]][comp[j]] += graph[i][j]
    return new_g, new_sz


# %%
def maek_U_V_Graph(U, V, sz):
    mark = np.zeros(sz, dtype=int)
    x = 1
    for i in range(sz):
        flag = 0
        if mark[i]:
            continue
        for j in range(sz):
            if U[i][j] == 1:
                flag += 1
        if flag > 0:
            dfs(U, mark, x, i, sz)
            x += 1
    new_graph, new_sz = make_components(V, mark, sz)
    return new_graph, new_sz


# %%
def laplacian(gamma, sz):
    L = np.zeros((sz - 1, sz - 1), dtype=int)
    for i in range(1, sz):
        for j in range(i, sz):
            if (i, j) in gamma.keys():
                L[i - 1][j - 1] = L[j - 1][i - 1] = -exp(gamma[(i, j)])
            if i == j:
                s = 0
                for k in range(sz):
                    if k == i:
                        continue
                    if (k, i) in gamma.keys():
                        s += exp(gamma[(k, i)])
                    if (i, k) in gamma.keys():
                        s += exp(gamma[(i, k)])
                L[i - 1][j - 1] = s
    return np.linalg.det(L)


# %%
def contract_Edge(graph, nodes, sz):
    new_g = np.zeros((sz - 1, sz - 1), dtype=int)
    for i in range(sz):
        for j in range(i + 1, sz):
            s = i
            t = j
            if nodes[0] < s:
                s -= 1
            if nodes[1] < s:
                s -= 1
            if nodes[0] < t:
                t -= 1
            if nodes[0] < t:
                t -= 1
            if (not (i in nodes)) and (not (j in nodes)):
                new_g[s][t] = graph[i][j]
                new_g[t][s] = graph[i][j]
            if (i == nodes[0]) and (j != nodes[1]):
                new_g[s][sz - 1] += graph[i][j]
                new_g[sz - 1][s] += graph[i][j]
            if (i != nodes[0]) and (j == nodes[1]):
                new_g[t][sz - 1] += graph[i][j]
                new_g[sz - 1][t] += graph[i][j]
    return new_g


# %%
def q(gamma, sz):
    return (laplacian(gamma, sz)) / (laplacian(gamma, sz))


# %%
def spanning_Tree_Distribution(z, sz):
    gamma = {}
    for i in range(sz):
        for j in range(i + 1, sz):
            gamma[(i, j)] = 0
    epsilon = 0.2
    while True:
        flag = 0
        for i in range(sz):
            for j in range(i + 1, sz):
                e = (i, j)
                q_e = q(gamma, sz)
                z_e = z[e]
                if q_e > (1 + epsilon) * z_e:
                    flag += 1
                    delta = ln((q_e * (1 - (1 + epsilon / 2) * z_e)) / ((1 - q_e) * (1 + epsilon / 2) * z_e))
                    gamma[e] -= delta
        if flag == 0:
            break
    return gamma


# %%
def sample_Tree(graph, gamma, sz):
    u = np.zeros((sz, sz), dtype=int)
    v = graph.copy()
    edges = []
    for i in range(sz):
        for j in range(i + 1, sz):
            if (i, j) in gamma.keys():
                edges.append(sz * i + j)
    random.shuffle(edges)
    cnt = 0
    for i in range(len(edges)):
        e = edges[i]
        g, new_sz = maek_U_V_Graph(u, v, sz)
        a = laplacian(gamma, new_sz)
        u[e[0]][e[1]] = 1
        g, new_sz = maek_U_V_Graph(u, v, sz)
        a_prime = laplacian(gamma, new_sz)
        u[e[0]][e[1]] = 0
        z = random.uniform(0, 1)
        if z <= (gamma[e] * a_prime / a):
            u[e[0]][e[1]] = 1
            cnt += 1
            if cnt == sz - 1:
                break
        else:
            v[e[0]][e[1]] = v[e[1]][e[0]] = 0
    return u


# %%
def direct_Tree(c, graph, sz):
    s = 0
    for i in range(sz):
        for j in range(i + 1, sz):
            u = i * sz + j
            v = j * sz + i
            if c[u] < c[v]:
                graph[i][j] = c[u]
                graph[j][i] = 0
                s += c[u]
            else:
                graph[j][i] = c[v]
                graph[i][j] = 0
                s += c[v]
    return graph, s


# %%
def make_demands(graph, sz):
    demands = []
    for i in range(sz):
        out = 0
        for j in range(sz):
            out += graph[i][j]
        inn = 0
        for j in range(sz):
            inn += graph[j][i]
        demands.append(out - inn)
    return demands


# %%
def make_Answer_From_X(c, x, graph, sz):
    z = {}
    for i in range(sz):
        for j in range(i + 1, sz):
            z[(i, j)] = (sz - 1) / sz * (x[i][j] + x[j][i])
    gamma = spanning_Tree_Distribution(z, sz)
    v = np.zeros((sz, sz), dtype=int)
    for i in range(sz):
        for j in range(i + 1, sz):
            if (i, j) in z.keys():
                v[i][j] = v[j][i] = 1
    cost = INF
    final = np.zeros((sz, sz), dtype=int)
    for i in range(math.ceil(6 * math.log(sz))):
        temp = sample_Tree(v, gamma, sz)
        temp_g, temp_cost = direct_Tree(c, temp, sz)
        if temp_cost < cost:
            cost = temp_cost
            final = temp_g
    demands = make_demands(final, sz)
    g = nx.DiGraph()
    for i in range(sz):
        g.add_node(i, demand=demands[i])
    for i in range(sz):
        for j in range(sz):
            g.add_edge(i, j, weight=graph[i][j])
    flow_dict = nx.min_cost_flow(g)
    return flow_dict


def make_initial_constraints(M, n):
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
    A_eq, b_eq, A_ub, b_ub, c, sz = make_initial_constraints(M, n)
    _, x = held_karp(A_eq, b_eq, A_ub, b_ub, c, sz)
    return make_Answer_From_X(c, np.array(x).reshape((sz, sz)), M, sz)
