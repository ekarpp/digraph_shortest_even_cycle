#!/usr/bin/env python3

def graph(p):
    grph = ""
    for i in range(p):
        elems = [
            (i + 1) % p,
            (i - 1 + p) % p,
            1
        ]
        for _ in range(p-2):
            elems[2] = (elems[2] * i) % p
        graph += " ".join(elems) + "\n"

import sys

for arg in sys.argv[1:]:
    graph(int(arg))
