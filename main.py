import sys

import FGM
from asadpour import asadpour
from utils import get_adjacency_matrix
f = open('nul', 'w')
sys.stderr = f
available_tests = ['br17.atsp', 'ft53.atsp', 'ft70.atsp', 'ftv170.atsp', 'ftv33.atsp', 'ftv35.atsp', 'ftv38.atsp',
                   'ftv44.atsp', 'ftv47.atsp', 'ftv55.atsp', 'ftv64.atsp', 'ftv70.atsp', 'kro124p.atsp', 'p43.atsp',
                   'rbg323.atsp', 'rbg358.atsp', 'rbg403.atsp', 'rbg443.atsp', 'ry48p.atsp']
solutions = dict(
    {"br17": 39, "ft53": 6905, "ft70": 38673, "ftv33": 1286, "ftv35": 1473, "ftv38": 1530, "ftv44": 1613, "ftv47": 1776,
     "ftv55": 1608, "ftv64": 1839, "ftv70": 1950, "ftv90": 1579, "ftv100": 1788, "ftv110": 1958, "ftv120": 2166,
     "ftv130": 2307, "ftv140": 2420, "ftv150": 2611, "ftv160": 2683, "ftv170": 2755, "kro124": 36230, "p43": 5620,
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
        print("FGM:")
        print(FGM.solve(M, n))
        print("Held-Karp lower bound")
        print(asadpour(M, n)[0])

    print("please enter the next test number")
    inp = input().strip()
