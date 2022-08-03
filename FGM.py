# %%
import numpy as np
from scipy import optimize

from constants import INF


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
    return nodes, comp


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
    new_sz = 0
    for i in range(sz):
        new_sz = max(new_sz, comp[i])
    new_sz += 1
    new_g = np.zeros((new_sz, new_sz), dtype=int)
    for i in range(new_sz):
        new_g[i][i] = INF
    in_edges = np.zeros(new_sz, dtype=int)
    for i in range(sz):
        for j in range(sz):
            if comp[i] != comp[j]:
                if graph[i][j] > new_g[comp[i]][comp[j]]:
                    new_g[comp[i]][comp[j]] = graph[i][j]
                    in_edges[comp[j]] = j
    return new_g, in_edges, new_sz


# %%
def solve(graph, sz):
    tsp_list = []
    nodes, comp = minimum_Matching(split_Graph(graph, sz), sz)
    new_graph, in_edges, new_sz = contract(graph, comp, sz)
    if new_sz == 1:
        return nodes, 0
    new_list, _ = solve(new_graph, new_sz)
    mark = np.zeros(sz, dtype=int)
    for i in range(new_sz):
        x = 0
        for j in range(sz):
            if comp[nodes[j]] == new_list[i]:
                x = j
                break
        while comp[nodes[x]] == new_list[i]:
            mark[nodes[x]] = 1
            tsp_list.append(nodes[x])
            x += 1
            if x == sz:
                break
    answer = 0
    for i in range(sz):
        if i == sz - 1:
            answer += graph[tsp_list[sz - 1]][tsp_list[0]]
        else:
            answer += graph[tsp_list[i]][tsp_list[i + 1]]
    return tsp_list, answer
