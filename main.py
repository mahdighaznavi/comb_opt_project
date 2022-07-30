import FGM
from asadpour import asadpour
from utils import get_adjacency_matrix

available_tests = ['br17.atsp']
print("Available tests:")
for i in range(len(available_tests)):
    print(str(i + 1) + ". " + available_tests[i])
inp = input().strip()
while inp != "exit":
    if inp.isdigit() and 1 <= inp.isdigit() <= len(available_tests):
        M, n = get_adjacency_matrix(available_tests[int(inp) - 1])
        print("FGM:")
        print(FGM.solve(M, n))
        print("Held-Karp lower bound")
        print(asadpour(M, n)[0])
        print()
        print("Available tests:")
        for i in range(len(available_tests)):
            print(str(i + 1) + ". " + available_tests[i])

    inp = input().strip()
