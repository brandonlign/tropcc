#!/usr/bin/env python3
"""Independent check that the mass-preserving body swaps g1=(1 2), g2=(3 4) are symmetries of the input:
 (i) the polynomial SET sys.json (c=1 not needed; exact over Z[a,b,c]) is mapped to itself up to nonzero rational scalars,
 (ii) the fan (rays, cones) is mapped to itself.
Then a certificate for one cell of an orbit transports to all cells of the orbit (apply the permutation to the identity).
Writes orbits.json: cell -> orbit representative (min cell index of the orbit)."""
import json,os,sys,itertools
from fractions import Fraction
import verify as V
rv,F=V.load_sys()   # F: list of {r-exponent: {(ea,eb): coeff}}  (c=1 substituted; AC homogeneous deg 1 so no info lost)
pairs=[(int(v[1]),int(v[2])) for v in rv]
def coordperm(g):   # g: body permutation dict; returns list P with new_exp[P[k]] = old_exp[k]
    idx={p:k for k,p in enumerate(pairs)}
    return [idx[tuple(sorted((g[i],g[j])))] for (i,j) in pairs]
def apply(P,e):
    out=[0]*len(e)
    for k,x in enumerate(e): out[P[k]]=x
    return tuple(out)
def normal(f):
    # scale so that the lexicographically first term's first coefficient is 1
    m0=min(f); ab0=min(f[m0]); s=Fraction(f[m0][ab0])
    return frozenset((m,ab,Fraction(c)/s) for m,d in f.items() for ab,c in d.items())
ok=True
Fn=set(normal(f) for f in F)
gens={'(12)':{1:2,2:1,3:3,4:4,5:5},'(34)':{1:1,2:2,3:4,4:3,5:5}}
for name,g in gens.items():
    P=coordperm(g)
    img=set(normal({apply(P,m):d for m,d in f.items()}) for f in F)
    print('polys invariant under',name,img==Fn); ok&=img==Fn
rays,cells=V.read_fan(os.path.join(V.N5,'fan.out'))
rix={tuple(r):i for i,r in enumerate(rays)}; cix={tuple(sorted(c)):k for k,c in enumerate(cells)}
maps=[]
for name,g in gens.items():
    P=coordperm(g)
    rm=[rix.get(apply(P,r)) for r in rays]
    good=None not in rm and all(tuple(sorted(rm[i] for i in c)) in cix for c in cells)
    print('fan invariant under',name,good); ok&=good
    maps.append([cix[tuple(sorted(rm[i] for i in c))] for c in cells] if good else None)
if ok:
    rep={}
    for k in range(len(cells)):
        if k in rep: continue
        orb={k}; st=[k]
        while st:
            x=st.pop()
            for mp in maps:
                if mp[x] not in orb: orb.add(mp[x]); st.append(mp[x])
        for y in orb: rep[y]=min(orb)
    json.dump(rep,open(os.path.join(V.HERE,'orbits.json'),'w'))
    print('cells',len(cells),'orbits',len(set(rep.values())))
sys.exit(0 if ok else 1)
