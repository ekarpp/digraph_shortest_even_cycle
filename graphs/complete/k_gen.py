#!/usr/bin/env python3

def K(n):
	l = [x for x in range(n)]
	grph = ""
	for i in range(n):
		lfilt = [str(x) for x in l if x != i]
		grph += f"{' '.join(lfilt)}\n"
	with open(f"k{n}", "w") as f:
		f.write(grph)

import sys

for arg in sys.argv[1:]:
	K(int(arg))
