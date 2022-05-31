#!/usr/bin/env python3

def exp(m):
    grph = ""
    for x in range(m):
        for y in range(m):
            nbors = [
                ((x+y)%m, y),
                (x, (y+x)%m),

                ((x-y+m)%m, y),
                (x, (y-x+m)%m),

                ((x+y+1)%m, y),
                (x, (y+x+1)%m),

                ((x-y+1+m)%m, y),
                (x, (x-y+1+m)%m)
            ]
            nbors = [str(xx + yy*m) for xx,yy in nbors]
            grph += " ".join(nbors) + "\n"
    with open(f"ee{m}", "w") as f:
        f.write(grph)

import sys

for arg in sys.argv[1:]:
    exp(int(arg))
