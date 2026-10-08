#!/usr/bin/env python3
"""Check that the gfan input data/supp.in has exactly the monomial supports of the 35 polynomials in sys.json
(with c=1: a,b,c treated as coefficients; a support monomial is present iff its coefficient polynomial in Q[a,b,c] is nonzero).
Pure-python parsing (no sympy for supp.in); sys.json via sympy Poly."""
import json,re,sys,os
import sympy as sp
H=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data')
D=json.load(open(os.path.join(H,'sys.json'))); rv=D['rvars']; R=sp.symbols(rv)
S1=[frozenset(sp.Poly(sp.sympify(f),*R).monoms()) for f in D['polys']]
txt=open(os.path.join(H,'supp.in')).read(); body=txt[txt.index('{')+1:txt.rindex('}')]
def mono(t):
    e=[0]*10
    for v,x in re.findall(r'(r\d\d)\^(\d+)',t): e[rv.index(v)]+=int(x)
    return tuple(e)
S2=[frozenset(mono(t) for t in p.split('+')) for p in body.split(',')]
ok=len(S1)==len(S2) and all(a==b for a,b in zip(S1,S2))
print('polys',len(S1),len(S2),'supports equal' if ok else 'MISMATCH'); sys.exit(0 if ok else 1)
