import sys

import FGM

import STV
from asadpour import asadpour
from utils import get_adjacency_matrix
from constants import Epsilon
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pylab

# f = open('nul', 'w')
# sys.stderr = f
available_tests = ['br17.atsp', 'ft53.atsp', 'ft70.atsp', 'ftv170.atsp', 'ftv33.atsp', 'ftv35.atsp', 'ftv38.atsp',
                   'ftv44.atsp', 'ftv47.atsp', 'ftv55.atsp', 'ftv64.atsp', 'ftv70.atsp', 'kro124p.atsp', 'p43.atsp',
                   'rbg323.atsp', 'rbg358.atsp', 'rbg403.atsp', 'rbg443.atsp', 'ry48p.atsp']
solutions = dict(
    {"br17": 39, "ft53": 6905, "ft70": 38673, "ftv33": 1286, "ftv35": 1473, "ftv38": 1530, "ftv44": 1613, "ftv47": 1776,
     "ftv55": 1608, "ftv64": 1839, "ftv70": 1950, "ftv90": 1579, "ftv100": 1788, "ftv110": 1958, "ftv120": 2166,
     "ftv130": 2307, "ftv140": 2420, "ftv150": 2611, "ftv160": 2683, "ftv170": 2755, "kro124p": 36230, "p43": 5620,
     "rbg323": 1326, "rbg358": 1163, "rbg403": 2465, "rbg443": 2720, "ry48p": 14422})
print(len(available_tests))
print(len(solutions))
print("Available tests:")
for i in range(len(available_tests)):
    print(str(i + 1) + ". " + available_tests[i])
inp = input().strip()
while inp != "exit":
    if inp.isdigit() and 1 <= inp.isdigit() <= len(available_tests):
        M, n = get_adjacency_matrix("tests/" + available_tests[int(inp) - 1])
        print("Solution: " + str(solutions[available_tests[int(inp) - 1][:-5]]))
        # STV.solve(M, n, Epsilon)
        # exit()
        print("FGM upper bound")
        vertex_list, answer = FGM.solve(M, n)
        print(answer)
        print(vertex_list)
        G = nx.DiGraph()
        nodes = np.arange(0, n).tolist()
        G.add_nodes_from(nodes)
        for i in range(len(vertex_list)):
            u = vertex_list[i]
            v = vertex_list[(i + 1) % len(vertex_list)]
            G.add_edges_from([(u, v)], weight=M[u][v])
        pos = nx.circular_layout(G)
        edge_labels = dict([((u, v,), d['weight'])
                            for u, v, d in G.edges(data=True)])
        nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), node_size=5)
        nx.draw_networkx_labels(G, pos, font_size=5)
        nx.draw_networkx_edges(G, pos)
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=5)
        pylab.show()
        print("Held-Karp lower bound:")
        tmp = asadpour(M, n)
        print("Asadpour upper bound")
        print(tmp)

    print("please enter the next test number")
    inp = input().strip()
