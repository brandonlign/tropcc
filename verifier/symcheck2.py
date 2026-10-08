#!/usr/bin/env python3
"""S17: independent exact check of the extra symmetry sigma = body permutation (13)(24) combined with the mass swap a<->b
(masses (a,a,b,b,1) -> (b,b,a,a,1) = the same multiset after relabelling bodies). Checks
 (i)  sys.json polynomial SET maps to itself up to nonzero rational scalars under r_ij -> r_sigma(i)sigma(j), (a,b)->(b,a);
 (ii) the augmented systems F+G and F+G+H (verify.augment/augmentH, U fixed) likewise;
 (iii) the fan (rays, cones) maps to itself.
Then for any certificate of cell k over Q(a,b) (an identity sum T_i g_i = L), applying sigma gives an identity for cell sigma(k)
with generators sigma(g_i) = (rational) * g_{pi(i)} and L(b,a) != 0, so cell sigma(k) is excluded for all (a,b) with L(b,a) != 0.
Since sigma (12) sigma = (34), the group <(12),(34),sigma> (dihedral, order 8) has orbits orb(k) u orb(sigma k).
Writes sigma.json: list, cell k -> sigma(k)."""
import json,os,sys
from fractions import Fraction
import verify as V
g={1:3,2:4,3:1,4:2,5:5}
def coordperm(rv):
    pairs=[(int(v[1]),int(v[2])) for v in rv[:10]]; idx={p:k for k,p in enumerate(pairs)}
    return [idx[tuple(sorted((g[i],g[j])))] for (i,j) in pairs]+list(range(10,len(rv)))
def apply(P,e):
    out=[0]*len(e)
    for k,x in enumerate(e): out[P[k]]=x
    return tuple(out)
def normal(f):
    m0=min(f); ab0=min(f[m0]); s=Fraction(f[m0][ab0])
    return frozenset((m,ab,Fraction(c)/s) for m,d in f.items() for ab,c in d.items())
def img(P,f): return {apply(P,m):{(eb,ea):c for (ea,eb),c in d.items()} for m,d in f.items()}
ok=True
rv,F=V.load_sys()
for name,(r,FF) in [('F',(rv,F)),('F+G',V.augment(rv,F)),('F+G+H',V.augmentH(rv,F))]:
    P=coordperm(r); A=sorted(map(normal,FF),key=str); B=sorted((normal(img(P,f)) for f in FF),key=str)
    good=set(A)==set(B) and len(set(A))==len(FF); print('poly set invariant under sigma:',name,good); ok&=good
rays,cells=V.read_fan(os.path.join(V.N5,'fan.out'))
P=coordperm(rv)
rix={tuple(x):i for i,x in enumerate(rays)}; cix={tuple(sorted(c)):k for k,c in enumerate(cells)}
rm=[rix.get(apply(P,x)) for x in rays]
good=None not in rm and all(tuple(sorted(rm[i] for i in c)) in cix for c in cells); print('fan invariant under sigma:',good); ok&=good
if ok:
    sm=[cix[tuple(sorted(rm[i] for i in c))] for c in cells]; assert all(sm[sm[k]]==k for k in range(len(cells)))
    json.dump(sm,open(os.path.join(V.HERE,'sigma.json'),'w')); print('wrote sigma.json, cells',len(sm))
sys.exit(0 if ok else 1)
