def get_adjacency_matrix(file_path):
    f = open(file_path, "r")
    f.readline()
    f.readline()
    f.readline()
    dimension = int(f.readline().split()[1])
    print(dimension)
    while f.readline().strip() != "EDGE_WEIGHT_SECTION":
        continue
    e = []
    for i in range(dimension):
        e.append([])
        while len(e[-1]) < dimension:
            e[-1].extend(map(int, f.readline().split()))

    f.close()
    return e, dimension
