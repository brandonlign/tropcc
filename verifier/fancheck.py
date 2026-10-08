#!/usr/bin/env python3
"""Cross-version fan comparison: compare the stored prevariety fan (data/fan.out, patched gfan 0.7 macOS)
with a separate cross-version recomputation (Debian gfan 0.6.2, aarch64) of the same input supp.in.
Compares primitive ray sets, and the sets of cones (as ray sets), dimension by dimension.
usage: fancheck.py OTHER_FAN.out"""
import sys,os
from math import gcd
from functools import reduce
import verify as V
def prim(r):
    g=reduce(gcd,[abs(x) for x in r]) or 1; return tuple(x//g for x in r)
def load(fn):
    rays,cells=V.read_fan(fn)
    R=[prim(r) for r in rays]
    return set(R), set(frozenset(R[i] for i in c) for c in cells)
A=load(os.path.join(V.N5,'fan.out')); B=load(sys.argv[1])
print('rays: stored',len(A[0]),'other',len(B[0]),'common',len(A[0]&B[0]))
print('cones: stored',len(A[1]),'other',len(B[1]),'common',len(A[1]&B[1]))
same=A[0]==B[0] and A[1]==B[1]; print('IDENTICAL' if same else 'DIFFERENT'); sys.exit(0 if same else 1)
