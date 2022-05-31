#!/usr/bin/env python3

def C(n):
    l = [str((x+1)%n) for x in range(n)]
    with open(f"c{n}", "w") as f:
        f.write(
            "\n".join(l) + "\n"
        )

import sys

for arg in sys.argv[1:]:
    C(int(arg))
