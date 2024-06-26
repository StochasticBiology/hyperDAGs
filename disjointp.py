import sys
from itertools import combinations as choose
from collections import defaultdict as D
from collections import namedtuple as T

import networkx as nx

Node = T("Node", "string e0 e1 dec")

from datetime import datetime as dt

INF = 2**31 - 1

"""
1. Read requirements
2. Add ZERO → s for marked s
3. Add edges
4. For requirement: if satisfied then remove
5. Transitive sparsification: if a→x and a→b and exists b→x path in REQ, remove a→x
6. For a→b in req: find free sink(a)→source(b) path
7. For a→b in req: find free a→b path
8. For a→b in req: find a→b path
"""


def node(x):
    e1 = sum(1 for e in x if e == "1")
    return Node(x, len(x) - e1, e1, int(x, 2))


def dist(s, t):
    u, v = s.string, t.string
    d = 0
    for i in range(len(u)):
        if u[i] == v[i]:
            continue
        if u[i] == "1":
            return float("inf")
        d += 1
    return d


def incomp(u, v):
    return dist(u, v) >= INF and dist(v, u) >= INF


def eb_val(k):
    return max(k - 1, 0)


def cost(H):
    I, O = degrees(H)
    V = set(s for s, _ in H.edges())
    return sum(eb_val(O[v]) for v in V)


def degrees(H):
    nodes = set()
    I = D(int)
    O = D(int)
    for s, t in H.edges():
        O[s] += 1
        I[t] += 1
    return I, O


def reach(G, s):
    for v in G.neighbors(s):
        yield v
        yield from reach(G, v)


def implicit_neighbor(s):
    string = s.string
    for i, e in enumerate(string):
        if e == "0":
            yield node(string[:i] + "1" + string[i + 1 :])


def implicit_path(s, t, H=None):
    # if H is not None we will try to find a path without branching
    if not dist(s, t) < 2**31:
        raise Exception()
    if s == t:
        return []

    for v in implicit_neighbor(s):
        if H and v in H.nodes() and not H.neighbors(v):
            if dist(v, t) < INF:
                p = implicit_path(v, t, H=H)
                if p is not None:
                    return [v] + p
        elif H and v in H.nodes():
            possible = H.neighbors(v)
            for u in possible:
                v = u
                if dist(v, t) < INF:
                    p = implicit_path(v, t, H=H)
                    if p is not None:
                        return [v] + p
        else:
            if dist(v, t) < INF:
                p = implicit_path(v, t)
                if p is not None:
                    return [v] + p


def sinks(N, s):
    """Given a graph N and a vertex s, output all sinks reachable from s"""
    if s not in N.nodes():
        return
    if not N.neighbors(s):
        yield s
    for v in N.neighbors(s):
        yield from sinks(N, v)


def clean(E, G):
    """Remove from E the requirements satisfied in G"""
    sat = set()
    V = set(G.nodes())
    for s, t in E:
        if s in V and t in V and nx.has_path(G, s, t):
            sat.add((s, t))
    retval = bool(sat)
    for r in sat:
        # print("Remove", r)
        E.remove(r)
    return retval


def main():
    start = dt.now()
    G = []
    for line in sys.stdin:
        e = line.split()
        a, b = node(e[0]), node(e[1])
        if dist(a, b) >= INF:
            a, b = b, a
        assert dist(a, b) < len(a.string), f"dist({a},{b}) = {dist(a,b)}"
        G.append((a, b))

    ZERO = node("0" * len(G[0][0].string))
    for s, t in list(G):
        G.append((ZERO, s))
        G.append((ZERO, t))

    E = set(G)
    H = nx.DiGraph()

    # ADD ALL EDGES
    for s, t in E:
        d = dist(s, t)
        if d == 1:
            H.add_edge(s, t)

    clean(E, H)

    update = True
    while update:
        update = False

        for s, t in E:
            ss = sinks(H, s)
            for x in ss:
                if dist(x, t) >= INF:
                    continue
                P = implicit_path(x, t, H)
                if P:
                    p = [x] + P
                    for i in range(len(P) - 1):
                        H.add_edge(P[i], P[i + 1])
                else:
                    P = implicit_path(x, t)
                    if P:
                        p = [x] + P
                        for i in range(len(P) - 1):
                            H.add_edge(P[i], P[i + 1])
        update = clean(E, H)

    for s, t in E:
        P = implicit_path(s, t)
        if not P:
            continue

        for i in range(len(P) - 1):
            H.add_edge(P[i], P[i + 1])

    clean(E, H)

    for s, t in E:
        P = implicit_path(s, t)
        if not P:
            continue

        for i in range(len(P) - 1):
            H.add_edge(P[i], P[i + 1])

    clean(E, H)

    print("Cost", cost(H))

    end = dt.now()
    ds = (end - start).total_seconds()
    print(round(ds * 1000, 3), "ms")

    with open("_out", "w") as fout:
        for s, t in H.edges():
            print(s.string, t.string, file=fout)


main()
