import csv
from itertools import combinations as choose
from collections import defaultdict as D
from collections import namedtuple as T

Node = T("Node", "string e0 e1 dec name L")
R = T("R", "source target")

ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
ALPHABET += ALPHABET.lower()
ALPHABET += "0123456798"
ALPHABET += "?!&#%-=+*^~'"
ALPHABET *= 10

from datetime import datetime as dt


def x_coords(V):
    y = D(list)
    for v in V:
        y[v.e1].append(v)
    for k in y:
        y[k] = sorted(y[k])
    x = {}
    for k in y:
        for i in range(len(y[k])):
            x[y[k][i]] = (i, len(y[k]))
    return x


def digraph_cost(G):
    return sum(eb_val(G.out_degree(v)) for v in G.nodes())


def node(x):
    name = "".join(ALPHABET[i] if x[i] == "1" else "_" for i in range(len(x)))
    e1 = sum(1 for e in x if e == "1")
    return Node(x, len(x) - e1, e1, int(x, 2), name, len(x))


def dist(s, t):
    u, v = s.string, t.string
    d = 0
    for i in range(len(u)):
        if u[i] == v[i]:
            continue
        if u[i] in (1, "1"):
            return float("inf")
        d += 1
    return d


def incomp(u, v):
    return dist(u, v) > 100 and dist(v, u) > 100


def eb_val(k):
    return max(k - 1, 0)


def cost(H):
    I, O = degrees(H)
    V = set(s for s, _ in H)
    return sum(eb_val(O[v]) for v in V)


def degrees(H):
    nodes = set()
    I = D(int)
    O = D(int)
    for s, t in H:
        O[s] += 1
        I[t] += 1
    return I, O


def reach(N, s):
    for v in N[s]:
        yield v
        yield from reach(N, v)


def path(N, s, t, deg=None):
    if deg is None:
        deg = D(list)
    L = len(s.string)
    if not dist(s, t) < L + 1:
        return None
    if s == t:
        return []
    for v in N[s]:
        if len(deg[v]) > 2:
            continue
        if dist(v, t) < L + 1:
            p = path(N, v, t)
            if p is not None:
                return [v] + p


def _OUT(v):
    s = v.string
    for i in range(v.L):
        if s[i] == "0":
            yield node(s[:i] + "1" + s[i + 1 :])


def implicit_path(L, s, t):
    assert L == len(s.string) == len(t.string)  # , f"{L=}, {s=}, {t=}"
    if s == t:
        return []
    if dist(s, t) > L:
        return None
    for v in _OUT(s):
        if dist(v, t) < L + 1:
            p = implicit_path(L, v, t)
            if p is not None:
                return [v] + p


def sinks(N, s):
    """Given a graph N and a vertex s, output all sinks reachable from s"""
    if not N[s]:
        yield s
    for v in N[s]:
        yield from sinks(N, v)


def clean_requirements(N, reqs):
    print("cleaning...")
    to_delete = set()
    for s, t in reqs:
        if t in reach(N, s):
            to_delete.add(R(s, t))
    for r in to_delete:
        reqs.remove(r)
    print("done")


def main(fname):
    start = dt.now()

    V = set()
    with open(fname, "r") as fin:
        rows = csv.reader(fin)
        G = []
        next(rows)
        for x, y in rows:
            s = node(x)
            t = node(y)
            V |= set((s, t))
            if dist(s, t) <= s.L:
                r = R(s, t)
            elif dist(t, s) <= s.L:
                r = R(t, s)
            else:
                print(
                    f"Impossible requirement, {s.string}, {t.string}", file=sys.stderr
                )
            G.append(r)
    L = s.L
    v_Z = node("0" * L)

    for v in V:
        if v.dec:
            G.append(R(v_Z, v))

    E = sorted(set(G))
    print(len(E), "requirements")
    H = set()
    for s, t in sorted(E):
        d = dist(s, t)
        if d == 1:
            H.add(R(s, t))

    for s, t in E:
        assert dist(s, t) <= L, f"{s.L}, {t.L}, {L}, {dist(s, t)}"

    I, O = degrees(H)
    nodes = set()
    for s, t in H:
        nodes.add(s)
        nodes.add(t)

    # Create neigh
    N = D(list)
    for s, t in H:
        N[s].append(t)

    clean_requirements(N, E)

    for s, t in E:
        ss = sinks(N, s)
        found = False
        for v in ss:
            if dist(v, t) < 100:
                found = True

    last_val = len(E) + 1
    while 1:
        if last_val == len(E):
            break
        last_val = len(E)

        for s, t in E:
            ss = sinks(N, s)
            for x in ss:
                P = implicit_path(L, x, t)
                if P:
                    p = [x] + P
                    for i in range(len(P) - 1):
                        H.add((P[i], P[i + 1]))
                        N[P[i]].append(P[i + 1])
                        N[P[i]] = sorted(set(N[P[i]]))
        clean_requirements(N, E)

    clean_requirements(N, E)

    for s, t in E:
        P = implicit_path(L, s, t)
        if not P:
            print("\n\n\nNO PATH!\n", L, s, t, dist(s, t))
        for i in range(len(P) - 1):
            H.add((P[i], P[i + 1]))
            N[P[i]].append(P[i + 1])
            N[P[i]] = sorted(set(N[P[i]]))

    clean_requirements(N, E)

    print("Cost", cost(H))

    end = dt.now()
    ds = (end - start).total_seconds()
    print(round(ds * 1000, 3), "ms")
    print(len(E), "remaining")
    for s, t in sorted(E):
        print("\tR:", s.string, t.string, dist(s, t))

    with open("g.csv", "w") as fout:
        print("source,target", file=fout)
        for s, t in H:
            print(f"{s.string},{t.string}", file=fout)

    return H


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        exit("usage: pypy3 disjointp.py data/mt-trans-manual.txt")
        exit("usage: pypy3 disjointp.py data/pt-trans-manual.txt")
    fname = sys.argv[1]
    H = main(fname)
