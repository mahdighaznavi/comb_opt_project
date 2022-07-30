def get_adjacency_matrix(file_path):
    f = open(file_path, "r")
    f.readline()
    f.readline()
    f.readline()
    dimension = int(f.readline().split()[1])
    while f.readline().strip() != "EDGE_WEIGHT_SECTION":
        continue
    e = []
    tmp = []
    for i in range(dimension):
        e.append([])
        for t in tmp:
            e[-1].append(t)
            if len(e[-1]) == dimension:
                break
        while len(e[-1]) < dimension:
            tmp = map(int, f.readline().split())
            for t in tmp:
                e[-1].append(t)
                if len(e[-1]) == dimension:
                    break
    f.close()
    return e, dimension
