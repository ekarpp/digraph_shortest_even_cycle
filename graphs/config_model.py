#!/usr/bin/env python3
import random

def create_adja(lst, n):
    adj = [[] for _ in range(n)]
    while lst:
        idx = random.randrange(0, len(lst))
        e1 = lst.pop(idx)
        idx = random.randrange(0, len(lst))
        e2 = lst.pop(idx)
        adj[e1].append(e2)
        adj[e2].append(e1)
    return adj

def adja_ok(adj):
    for i in range(len(adj)):
        if i in adj[i] or len(adj[i]) != len(set(adj[i])):
            return False
    return True

def cfg_model(n, d):
    W = []
    for i in range(n):
        for j in range(d):
            W.append(i)
    G = create_adja(W.copy(), n)
    while not adja_ok(G):
        G = create_adja(W.copy(), n)
    graph = "\n".join(
        [" ".join([str(x) for x in l]) for l in G]
    ) + "\n"
    with open(f"cm{n}_{d}", "w") as f:
        f.write(graph)

import sys

n = int(sys.argv[1])
d = int(sys.argv[2])

if n*d % 2 == 0:
    cfg_model(n, d)
