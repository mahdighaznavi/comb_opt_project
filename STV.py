import math

from held_karp import held_karp
import numpy as np
from scipy import optimize
from constants import INF


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


def epsilon_prime(epsilon, alpha):
    return epsilon / (3 + 4 * alpha + 1 / (2 * alpha))


def l(v, v_b, alpha, kappa, beta, epsilon, lp, L, y):
    if type(v) == list:
        ret = 0
        for vv in v:
            ret += l(vv, v_b, alpha, kappa, beta, epsilon, lp, L, y)
        return ret
    s = 0
    for u in range(L.shape[1]):
        if u not in v_b:
            for i in range(L.shape[0]):
                if sum(L[i]) == 1 and L[i][u] == 1:
                    s += y[i]
                    break
    if v in v_b:
        return (kappa * lp + beta * 2 * s) / len(v_b)
    else:
        y_v = 0
        for i in range(L.shape[0]):
            if sum(L[i]) == 1 and L[i][v] == 1:
                y_v = y[i]
        return (1 + epsilon_prime(epsilon, alpha)) * 2 * alpha * 2 * y_v + (
                epsilon_prime(epsilon, alpha) / L.shape[1]) * s


def dfs(g, v, mark, cnt, mask, v_b):
    if mask[v]:
        return mark
    mask[v] = True
    mark[v] = cnt
    for i in range(len(g[v])):
        if i not in v_b:
            mark = dfs(g, g[v][i][0], mark, cnt, mask, v_b)

    return mark


def get_components(g, v_b):
    mark = []
    for i in range(len(g)):
        mark.append(-1)

    cnt = 0
    for i in range(len(g)):
        if (i not in v_b) and mark[i] == -1:
            mask = np.zeros(len(g))
            mark = dfs(g, i, mark, cnt + 1, mask, v_b)
            cnt += 1

    ret = []

    t = list(set(mark))
    t.remove(-1)
    mark_num = dict({})
    for i in range(len(t)):
        mark_num[t[i]] = i
        ret.append([])

    for i in range(len(g)):
        if i not in v_b:
            ret[mark_num[mark[i]]].append(i)

    return ret


def union(G, H, B):
    new_H = []
    for i in range(len(H)):
        new_H.append([])
        for j in range(len(H[i])):
            new_H.append(H[i][j])
    for b in B:
        new_H[b[0]].append(G[b[0]][b[1]])

    return new_H


def is_connected(g):
    mark = [-1] * len(g)
    cnt = 0
    for i in range(len(g)):
        if mark[i] == -1:
            mask = np.zeros(len(g))
            mark = dfs(g, i, mark, cnt, mask, [])
            cnt += 1

    return len(list(set(mark))) == 1


def subtour_cover(I, B, H):
    G, L, x, y, active = I
    V_B = []
    for e in B:
        V_B.append(e[0])
    V_B = list(set(V_B))
    l = [np.ones(len(G))]
    L_size = []
    for i in range(len(L)):
        L_size.append((sum(L[i]), i))
    L_size.sort()
    for i in range(len(L_size) - 1, -1, -1):
        if L_size[i][0] > 1:
            l.append(L[L_size[i][1]])
    r = [-1] * len(G)
    for i in range(len(G)):
        for j in range(len(l) - 1, -1, -1):
            if l[j][i] == 1:
                r[i] = j
                break
    if min(r) < 0:
        print("Error in sub tour cover 1")
    num_edges = 0
    edge_num = []
    cnt = 0
    G_01 = []
    for i in range(2 * len(G)):
        G_01.append([])
    for i in range(len(G)):
        num_edges += len(G[i])
        edge_num.append([])
        G_01[i + len(G)].append((i, 0))
        if i in V_B:
            G_01[i].append((i + len(G), 0))
        for j in range(len(G[i])):
            edge_num[-1].append(cnt)
            cnt += 1
            if r[i] <= r[G[i][j][0]]:
                G_01[i].append((G[i][j][0], G[i][j][1]))
            if r[i] >= r[G[i][j][0]]:
                G_01[i + len(G)].append((G[i][j][0] + len(G), G[i][j][1]))

    bounds = []

    A_ub = []
    b_ub = []
    eq_num = []
    eq_cnt = 0
    for i in range(len(x)):
        if i not in V_B:
            A_ub.append(np.zeros(num_edges))
            b_ub.append(0)
            eq_num.append(eq_cnt)
            eq_cnt += 1
        else:
            eq_num.append(None)
    for i in range(len(x)):
        for j in range(len(x[i])):
            if i not in V_B:
                A_ub[eq_num[i]][edge_num[i][j]] = -1
            if x[i][j][0] not in V_B:
                A_ub[eq_num[x[i][j][0]]][edge_num[i][j]] = 1
            if r[i] < r[x[i][j][0]]:
                bounds.append((x[i][j][1], x[i][j][1]))
            elif r[i] > r[x[i][j][0]]:
                bounds.append((0, 0))
            else:
                bounds.append((0, x[i][j][1]))
    c = np.zeros(num_edges)
    W = get_components(H, V_B)
    component = -1 * np.ones(len(G))
    for i in range(len(W)):
        for v in W:
            component[v] = i
    for i in range(len(G)):
        for j in range(len(G[i])):
            if component[i] != component[G[i][j][0]]:
                c[edge_num[i][j]] += 1
    b_ub = np.array(b_ub)
    A_ub = np.array(A_ub)
    c = np.array(c)
    print(b_ub.shape)
    print(A_ub.shape)
    print(b_ub.shape)
    f_tilde = optimize.linprog(c=c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='interior-point', )

    c = np.ones(num_edges)
    for i in range(num_edges):
        bounds[i] = (bounds[i][0], f_tilde.x[i])

    f = optimize.linprog(c=c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='interior-point', )

    G_f = []


def has_intersect(a1, a2):
    for e in a1:
        if e in a2:
            return True
    return False


def phi(W_tildes, v_b, alpha, kappa, beta, epsilon, lp, L, y):
    ret = 0
    for w in W_tildes:
        ret += l(w, v_b, alpha, kappa, beta, epsilon, lp, L, y) ** (
                1 + math.log((2 + epsilon_prime(epsilon, alpha)) / epsilon_prime(epsilon, alpha),
                             base=1 + epsilon_prime(epsilon, alpha)))


def Svensson(I, B, H_tilde, alpha, kappa, beta, epsilon, lp):
    G, L, x, y, active = I
    W_tilde = []
    V_B = []
    for e in B:
        V_B.append(e[0])
    V_B = list(set(V_B))
    W_tilde_tmp = get_components(H_tilde, V_B)
    for i in range(len(W_tilde_tmp)):
        W_tilde_tmp[i] = (l(W_tilde_tmp[i], V_B, alpha, kappa, beta, epsilon, lp, L, y), W_tilde_tmp[i])
    W_tilde_tmp.sort()
    for i in range(len(W_tilde_tmp) - 1, -1, -1):
        W_tilde.append(W_tilde_tmp[i][1])

    H = H_tilde

    tmp = union(G, H, B)
    while not is_connected(tmp):
        F_prime = subtour_cover(I, B, H)
        F_prime_comps = get_components(F_prime, [])
        tmp_comps = get_components(tmp, [])
        for i in range(len(F_prime_comps)):
            remove_edges = False
            for j in range(len(tmp_comps)):
                if all(x in tmp_comps[j] for x in F_prime_comps[i]):
                    remove_edges = True
                    break
            if remove_edges:
                for j in range(len(F_prime_comps[i])):
                    F_prime[F_prime_comps[i][j]] = []
        F_comps = get_components(F_prime, [])
        F = []
        c_F = []
        for i in range(len(W_tilde)):
            F.append([])
            c_F.append(0)
        indices = []
        for comp in F_comps:
            indices.append(None)
            for j in range(len(W_tilde)):
                if has_intersect(comp, W_tilde[j]):
                    indices[-1] = j
                    for k in range(len(comp)):
                        for i in range(len(F_prime[comp[k]])):
                            F[j].append((comp[k], i))
                            # c_F[j] +=
                    break

        # for i in range(len(W_tilde)):
        #     for j in range()


def sum_subset_neq(v, W, L, y, active):
    ret = 0
    for k in range(len(L)):
        if active[k] and L[k][v] == 1:
            subset = True
            not_eq = False
            for p in range(len(W)):
                if W[p] != L[k][p]:
                    not_eq = True
                if W[p] == 0 and L[k][p] == 1:
                    subset = False
            if subset and not_eq:
                ret += y[k]
    return ret


def D(W, L, y, p, active):
    mx = -INF
    u_star = -1
    v_star = -1
    for i in range(len(W)):
        if W[i] == 1:
            for j in range(len(W)):
                if i != j and W[j] == 1:
                    s = 0
                    check = False
                    for k in range(len(p[i])):
                        if p[i][k][0] == j:
                            s = p[i][k][1]
                            check = True
                            break
                    if not check:
                        print("Bug in D")
                    s += sum_subset_neq(i, W, L, y, active)
                    s += sum_subset_neq(j, W, L, y, active)
                    if s > mx:
                        mx = s
                        u_star = i
                        v_star = j

    return mx, u_star, v_star


def reduce_to_vertebrate_pairs(I, W, epsilon, lp):
    G, L, x, y, active = I
    G_prime = []
    reverse_G_prime = []
    L_prime = []
    new_y = []
    new_x = []
    new_active = []
    ret = []

    DW, u_star, v_star = D(W, L, y, G, active)
    print(u_star)
    print(v_star)
    print(DW)
    print(L)
    print(active)
    contract = []
    G_prime_idx = [-1] * len(G)
    set_to_contract = []
    for i in range(len(L)):
        if active[i] and L[i][u_star] == False and L[i][v_star] == False:
            subset = True
            for j in range(len(W)):
                if W[j] == 0 and L[i][j] == 1:
                    subset = False
                    break

            if subset:
                set_to_contract.append((sum(L[i]), i))
    set_to_contract.sort()
    for i in range(len(set_to_contract) - 1, -1, -1):
        idx = set_to_contract[i][1]
        new = 0
        for j in range(len(L[idx])):
            if L[idx][j] == 1:
                if G_prime_idx[j] == -1:
                    G_prime_idx[j] = len(contract)
                    if new == -1:
                        print("Error in contraction 1")
                        exit()
                    new = 1
                else:
                    if new == 1:
                        print("Error in contraction 2")
                        exit()
                    new = -1
        if new == 1:
            ret.extend(reduce_to_vertebrate_pairs(I, L[idx], epsilon, lp))
            contract.append(L[idx])
            for j in range(len(L_prime)):
                L_prime[j].append(0)
            L_prime.append([])
            for j in range(len(L_prime)):
                L_prime[-1].append(0)
            L_prime[-1][-1] = 1
            new_y.append(y[idx] + D(L[idx], L, y, G, active)[0] / 2)
            new_active.append(True)
    if sum(W) != len(W):
        tmp = []
        for i in range(len(W)):
            tmp.append(0)
            if W[i] != 1:
                G_prime_idx[i] = len(contract)
                tmp[i] = 1
        contract.append(np.array(tmp))
        for i in range(len(L_prime)):
            L_prime[i].append(0)
    for i in range(len(W)):
        if G_prime_idx[i] == -1:
            G_prime_idx[i] = len(contract)
            contract.append(np.zeros(len(G)))
            contract[-1][i] = 1
            for j in range(len(L_prime)):
                L_prime[j].append(0)
    if sum(W) != len(W):
        L_prime.append([])
        for i in range(len(contract)):
            L_prime[-1].append(0)
        for i in range(len(W)):
            if W[i] == 1:
                L_prime[-1][G_prime_idx[i]] = 1
        new_y.append(DW / 2)
        new_active.append(True)

    for i in range(len(L)):
        if active[i] and (L[i][u_star] == 1 or L[i][v_star] == 1):
            subset = True
            not_eq = False
            for j in range(len(W)):
                if W[j] != L[i][j]:
                    not_eq = True
                if W[j] == 0 and L[i][j] == 1:
                    subset = False
                    break
            if subset and not_eq:
                L_prime.append([])
                for j in range(len(contract)):
                    L_prime[-1].append(0)
                for j in range(len(L[i])):
                    if L[i][j] == 1:
                        L_prime[-1][G_prime_idx[j]] = 1
                new_y.append(y[i])
                new_active.append(True)
    for i in range(len(contract)):
        G_prime.append([])
        reverse_G_prime.append([])
        new_x.append([])
        for v in range(len(contract[i])):
            if contract[i][v] == 1:
                for j in range(len(G[v])):
                    G_prime[-1].append((G_prime_idx[G[v][j][0]], G[v][j][1]))
                    reverse_G_prime[-1].append((v, j))
                    new_x[-1].append((G_prime_idx[G[v][j][0]], x[v][j][1]))

    B = []
    for i in range(len(G_prime[G_prime_idx[u_star]])):
        if G_prime[G_prime_idx[u_star]][i][0] == G_prime_idx[v_star]:
            B.append((G_prime_idx[u_star], i))
            break
    for i in range(len(G_prime[G_prime_idx[v_star]])):
        if G_prime[G_prime_idx[v_star]][i][0] == G_prime_idx[u_star]:
            B.append((G_prime_idx[v_star], i))
            break
    ret.extend(B)
    I_prime = (G_prime, np.array(L_prime), new_x, np.array(new_y), new_active)

    H_tilde = []
    for i in range(len(contract)):
        H_tilde.append([])

    F_prime = Svensson(I_prime, B, H_tilde, 3, 2, 1, epsilon, lp)
    for i in range(len(F_prime)):
        for j in range(len(F_prime[i])):
            ret.append(reverse_G_prime[i][F_prime[i][j]])

    return ret


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


def solve(M, n, epsilon):
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
    new_L, new_y, new_nodes, num_members, active = find_laminar_y(L, y, nodes)
    G = []
    new_x = []
    for i in range(n):
        G.append([])
        new_x.append([])
        for j in range(n):
            if i != j:
                G[-1].append((j, M[i][j]))
                new_x[-1].append((j, x_star.x[i * n + j]))

    F = reduce_to_vertebrate_pairs((G, np.array(new_nodes[1:]), new_x, new_y[1:], active[1:]), np.ones(n), epsilon,
                                   x_star.fun)
