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
    cnt = 0
    for i in range(sz):
        if comp[i] == 0:
            cnt += 1
        new_sz = max(new_sz, comp[i])
    new_sz += cnt
    ind = np.zeros(sz, dtype=int)
    bef = 0
    for i in range(sz):
        if comp[i] != 0:
            bef += 1
            ind[i] = cnt + comp[i] - 1
        else:
            ind[i] = i - bef
    new_g = np.zeros((new_sz, new_sz), dtype=int)
    for i in range(sz):
        for j in range(i + 1, sz):
            if not ((i, j) in graph.keys()):
                continue
            if comp[i] == 0 and comp[j] == 0:
                new_g[ind[i]][ind[j]] -= exp(graph[(i, j)])
                new_g[ind[i]][ind[i]] += exp(graph[(i, j)])
                new_g[ind[j]][ind[j]] += exp(graph[(i, j)])
            elif comp[i] != comp[j]:
                new_g[ind[i]][ind[j]] -= exp(graph[(i, j)])
                new_g[ind[i]][ind[i]] += exp(graph[(i, j)])
                new_g[ind[j]][ind[j]] += exp(graph[(i, j)])
    temp = np.zeros((new_sz - 1, new_sz - 1), dtype=int)
    for i in range(new_sz - 1):
        for j in range(new_sz - 1):
            temp[i][j] = new_g[i + 1][j + 1]
    return np.linalg.det(temp)


# %%
def make_U_V_Graph(gamma, U, V, sz):
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
    dic_graph = {}
    for i in range(sz):
        for j in range(i + 1, sz):
            if V[i][j] == 1:
                dic_graph[(i, j)] = gamma[(i, j)]
    new_graph = make_components(dic_graph, mark, sz)
    return new_graph


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
                    if (i, k) in gamma.keys():
                        s += exp(gamma[(i, k)])
                    if (k, i) in gamma.keys():
                        s += exp(gamma[(k, i)])
                L[i - 1][j - 1] = s
    return np.linalg.det(L)


# %%
def contract_Edge(graph, nodes, sz):
    new_graph = np.zeros((sz - 1, sz - 1), dtype=int)
    for i in range(sz):
        for j in range(i + 1, sz):
            if not ((i, j) in graph.keys()):
                continue
            if (i == nodes[0]) and (j == nodes[1]):
                continue
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
                new_graph[s][t] -= exp(graph[(i, j)])
                new_graph[s][s] += exp(graph[(i, j)])
                new_graph[t][t] += exp(graph[(i, j)])
            elif (i == nodes[0]) and (j != nodes[1]):
                new_graph[t][sz - 2] -= exp(graph[(i, j)])
                new_graph[t][t] += exp(graph[(i, j)])
                new_graph[sz - 2][sz - 2] += exp(graph[(i, j)])
            elif (i != nodes[0]) and (j == nodes[1]):
                new_graph[s][sz - 2] -= exp(graph[(i, j)])
                new_graph[s][s] += exp(graph[(i, j)])
                new_graph[sz - 2][sz - 2] += exp(graph[(i, j)])
            elif i == nodes[1]:
                new_graph[t][sz - 2] -= exp(graph[(i, j)])
                new_graph[t][t] += exp(graph[(i, j)])
                new_graph[sz - 2][sz - 2] += exp(graph[(i, j)])
            elif j == nodes[0]:
                new_graph[s][sz - 2] -= exp(graph[(i, j)])
                new_graph[s][s] += exp(graph[(i, j)])
                new_graph[sz - 2][sz - 2] += exp(graph[(i, j)])
    for i in range(sz - 1):
        for j in range(i + 1, sz - 1):
            new_graph[j][i] = new_graph[i][j]
    temp = np.zeros((sz - 2, sz - 2), dtype=int)
    for i in range(sz - 2):
        for j in range(sz - 2):
            temp[i][j] = new_graph[i + 1][j + 1]
    return np.linalg.det(temp)


# %%
def q(e, gamma, sz):
    return (contract_Edge(gamma, e, sz)) / (laplacian(gamma, sz))


# %%
def spanning_Tree_Distribution(z, sz):
    gamma = {}
    epsilon = 0.2
    for i in range(sz):
        for j in range(i + 1, sz):
            if not ((i, j) in z.keys()):
                continue
            gamma[(i, j)] = 0
    while True:
        flag = 0
        for i in range(sz):
            for j in range(i + 1, sz):
                e = (i, j)
                if not (e in gamma.keys()):
                    continue
                q_e = q(e, gamma, sz)
                z_e = z[e]
                if q_e > ((1 + epsilon) * z_e):
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
        a = make_U_V_Graph(gamma, u, v, sz)
        u[e // sz][e % sz] = 1
        a_prime = make_U_V_Graph(gamma, u, v, sz)
        u[e // sz][e % sz] = 0
        z = random.uniform(0, 1)
        if z <= (1 * a_prime / a):
            u[e // sz][e % sz] = 1
            cnt += 1
            if cnt == sz - 1:
                break
        else:
            v[e // sz][e % sz] = v[e % sz][e // sz] = 0
    return u


# %%
def direct_Tree(c, graph, sz):
    s = 0
    for i in range(sz):
        for j in range(i + 1, sz):
            if graph[i][j] == 0:
                continue
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
            if (x[i * sz + j] + x[j * sz + i]) < 1e-9:
                continue
            z[(i, j)] = (sz - 1) / sz * (x[i * sz + j] + x[j * sz + i])
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
        if len(temp) < sz - 1:
            continue
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
    flow_dict = nx.min_cost_flow_cost(g)
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
    return make_Answer_From_X(c, x, M, sz)
