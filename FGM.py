# %%
import numpy as np
from scipy import optimize
from constants import INF


def minimum_Matching(graph, sz):
    row_ind, col_ind = optimize.linear_sum_assignment(graph)
    nodes = []
    mark = np.zeros(sz, dtype=int)
    comp = np.zeros(sz, dtype=int)
    x = 0
    for i in range(sz):
        t = i
        if (mark[t] == 0):
            x += 1
        while (not mark[t]):
            comp[t] = x - 1
            nodes.append(t)
            mark[t] = 1
            t = col_ind[t] - sz
    cost = graph[row_ind, col_ind].sum()
    return cost, nodes, comp


# %%
def split_Graph(graph, sz):
    new_Graph = np.zeros((2 * sz, 2 * sz), dtype=int)
    for i in range(sz):
        for j in range(sz):
            new_Graph[i][j] = INF
    for i in range(sz):
        for j in range(sz):
            new_Graph[i][sz + j] = graph[i][j]
    return new_Graph


# %%
def contract(graph, comp, sz):
    newSz = 0
    for i in range(sz):
        newSz = max(newSz, comp[i])
    newSz += 1
    newG = np.zeros((newSz, newSz), dtype=int)
    for i in range(newSz):
        newG[i][i] = INF
    in_Edges = np.zeros(newSz, dtype=int)
    for i in range(sz):
        for j in range(sz):
            if comp[i] != comp[j]:
                if (graph[i][j] > newG[comp[i]][comp[j]]):
                    newG[comp[i]][comp[j]] = graph[i][j]
                    in_Edges[comp[j]] = j
    return newG, in_Edges, newSz


# %%
def solve(graph, sz):
    tsp_List = []
    cost, nodes, comp = minimum_Matching(split_Graph(graph, sz), sz)
    newGraph, in_Edges, newSz = contract(graph, comp, sz)
    if (newSz == 1):
        return nodes, cost
    added = 0
    mark = np.zeros(sz, dtype=int)
    for i in range(newSz):
        if i == newSz - 1:
            added += graph[in_Edges[i]][in_Edges[0]]
        else:
            added += graph[in_Edges[i]][in_Edges[i + 1]]
        x = in_Edges[i]
        st = 0
        en = 0
        for j in range(sz):
            if nodes[j] == x:
                st = j
                en = j
        while (comp[nodes[st]] == i):
            mark[nodes[st]] = 1
            tsp_List.append(nodes[st])
            st += 1
            if (st == sz):
                break
        for j in range(sz):
            if comp[nodes[j]] == i:
                en = j
                break
        while (mark[nodes[en]] == 0):
            mark[nodes[en]] = 1
            tsp_List.append(nodes[en])
            en += 1
    new_List, newCost = solve(newGraph, newSz)
    return tsp_List, cost + added + newCost
