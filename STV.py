from held_karp import held_karp
import numpy as np
from scipy import optimize


def make_initial_constraints(M, n):
    c = []
    for i in range(n):
        for j in range(n):
            c.append(M[i][j])

    b_eq = np.zeros(n, dtype=int)
    A_eq = np.zeros((n, n * n), dtype=int)

    for i in range(n):
        b_eq[i] = 0
        for j in range(n):
            A_eq[i][i * n + j] += 1
            A_eq[i][j * n + i] -= 1

    b_ub = []
    A_ub = []
    b_ub.append(0)
    A_ub.append([])
    for i in range(n * n):
        A_ub[-1].append(0)

    return A_eq, b_eq, A_ub, b_ub, c


def check_laminar(l1, l2):
    """

    :param l1: first column
    :param l2: second column
    :return: None if crossing, 1 if l1 in l2, -1 if l2 in l1, 0 otherwise
    """
    one_m_tow = False
    two_m_one = False
    intersect = False
    for i in range(len(l1)):
        if l1[i] == 1 and l2[i] == 1:
            intersect = True
        elif l1[i] == 0 and l2[i] == 1:
            two_m_one = True
        elif l1[i] == 1 and l2[i] == 0:
            one_m_tow = True
    if intersect and one_m_tow and two_m_one:
        return None
    elif one_m_tow and two_m_one:
        return 0
    elif one_m_tow:
        return -1
    elif two_m_one:
        return 1


def minus(l1, l2):
    ret = np.zeros_like(l1)
    for i in range(len(l1)):
        if l1[i] == 1 and l2[i] == 0:
            ret[i] = 1
    return ret


def binary_search(t_array, x, l, r):
    if l + 1 == r:
        return l
    mid = (l + r) // 2
    if t_array[mid][0] <= x:
        return binary_search(t_array, x, mid, r)
    else:
        return binary_search(t_array, x, l, mid)


def check_equal(l1, l2):
    for i in range(len(l1)):
        if l1[i] != l2[i]:
            return False

    return True


def find_laminar_y(L, y, nodes):
    pre_check = []
    # debug
    print("Deb")
    print(L[:, 17])
    print(nodes)
    for i in range(L.shape[1]):
        pre_check.append(0)
        for j in range(len(y)):
            pre_check[-1] += y[j] * L[j][i]
    active = []
    num_members = []
    base_members = []
    for i in range(len(nodes)):
        active.append(False)
        num_members.append((sum(nodes[i]), i))
    num_members.sort()
    for i in range(len(num_members)):
        base_members.append((num_members[i][0], num_members[i][1]))

    def add(idx):
        nonlocal L
        nonlocal y
        nonlocal nodes
        activate = True
        repeat = True
        while repeat:
            repeat = False
            for n in range(len(num_members) - 1, -1, -1):
                tmp = num_members[n][1]
                if active[tmp]:
                    check = check_laminar(nodes[tmp], nodes[idx])
                    if check is None:
                        if y[idx] <= y[tmp]:
                            m = y[idx]
                            y[idx] = 0.0
                            y[tmp] -= m
                            activate = False
                        else:
                            m = y[tmp]
                            y[tmp] = 0.0
                            y[idx] -= m
                            active[tmp] = False
                        new_m_old = minus(nodes[idx], nodes[tmp])
                        old_m_new = minus(nodes[tmp], nodes[idx])
                        a = []
                        if sum(new_m_old) <= sum(old_m_new):
                            a.append(new_m_old)
                            a.append(old_m_new)
                        else:
                            a.append(old_m_new)
                            a.append(new_m_old)
                        should_add = [None, None]
                        for t in range(2):
                            num_ones = sum(a[t])
                            row = binary_search(num_members, num_ones, 0, len(num_members))
                            exist = False
                            while row >= 0 and num_members[row][0] == num_ones:
                                if check_equal(a[t], nodes[row]):
                                    y[row] += m
                                    if not active[row]:
                                        should_add[t] = row
                                    exist = True
                                    break
                                row -= 1
                            if not exist:
                                should_add[t] = len(nodes)
                                num_members.append((num_ones, len(nodes)))
                                num_members.sort()
                                nodes = np.concatenate((nodes, np.reshape(a[t], (1, a[t].shape[0]))))
                                L = np.concatenate((L, np.zeros((1, L.shape[1]))))
                                for e in range(L.shape[1]):
                                    u = e // nodes.shape[1]
                                    v = e % nodes.shape[1]
                                    if a[t][u] == 1 and a[t][v] == 0:
                                        L[-1][e] = 1
                                y = np.concatenate((y, np.array([m])))
                                active.append(False)

                        for t in range(2):
                            if should_add[t] is not None:
                                repeat = True
                                add(should_add[t])

                if repeat or not activate:
                    break

        if activate:
            active[idx] = True

    for i in range(len(base_members)):
        if not active[base_members[i][1]] and y[base_members[i][1]] > 0.0:
            add(base_members[i][1])

    return L, y, nodes, num_members, active


def solve(M, n):
    A_eq, b_eq, A_ub, b_ub, c = make_initial_constraints(M, n)
    x_star, A_ub, b_ub, nodes = held_karp(A_eq, b_eq, A_ub, b_ub, c, n, True)
    sym_b_ub = 2 * np.array(b_ub)
    A_ub = np.array(A_ub)
    for i in range(A_ub.shape[0]):
        for j in range(A_ub.shape[1]):
            if A_ub[i][j] == -1:
                u = j // n
                v = j % n
                A_ub[i][v * n + u] = -1
    A_eq = np.array(A_eq)
    b_eq = np.array(b_eq)
    dual_weight = np.concatenate((b_eq, -b_eq, sym_b_ub))
    dual_matrix = np.concatenate((-A_eq.transpose(), A_eq.transpose(), -np.array(A_ub).transpose()), axis=1)
    dual_right_side = np.array(c)
    y_star = optimize.linprog(c=dual_weight, A_ub=dual_matrix, b_ub=dual_right_side, method='interior-point', )
    y = y_star.x[-len(A_ub):]
    a = y_star.x[:len(A_eq)] - y_star.x[len(A_eq):2 * len(A_eq)]
    L = -np.array(A_ub)
    for i in range(L.shape[0]):
        print(nodes[i])
        print(L[i])
    new_L, new_y, new_nodes, num_members, active = find_laminar_y(L, y, nodes)
