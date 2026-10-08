"""Exact Sturm sign-variation certificate for the degree-eight mass factor."""
from pathlib import Path
import json,sympy as s
ROOT=Path(__file__).resolve().parents[1];x=s.symbols('x');P=x**8-x**7+x**6-x**4+8*x**3+3*x**2-3*x+1
seq=[s.Poly(P,x,domain=s.QQ),s.Poly(s.diff(P,x),x,domain=s.QQ)]
while seq[-1].degree()>0:seq.append(-seq[-2].rem(seq[-1]))
# Verify every Euclidean remainder relation explicitly.
for i in range(2,len(seq)):
 q,r=seq[i-2].div(seq[i-1]);assert seq[i-2]==q*seq[i-1]-seq[i]
signs0=[s.sign(f.eval(0)) for f in seq];signsinf=[s.sign(f.LC()) for f in seq]
def variations(values):
 z=[v for v in values if v];return sum(u!=v for u,v in zip(z,z[1:]))
assert variations(signs0)-variations(signsinf)==0;assert P.subs(x,0)==1
out={'polynomial':str(P),'sturm_sequence':[str(f.as_expr()) for f in seq],'signs_at_zero':list(map(int,signs0)),'signs_at_positive_infinity':list(map(int,signsinf)),'variation_counts':[variations(signs0),variations(signsinf)],'positive_roots':0,'conclusion':'P(x)>0 for x>0. For positive a,c, the mass resultant for cell 346 vanishes only at a=4c.'}
(ROOT/'results/sturm_mass_factor.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='sturm_sequence'})
