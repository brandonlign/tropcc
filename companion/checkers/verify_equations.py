"""Independently derive published distance equations, compare up to Laurent units.
Formula source: Jensen-Leykin arXiv:2301.02305v2, section 2.2.
"""
from pathlib import Path
import itertools,json,time
import sympy as s
ROOT=Path(__file__).resolve().parents[1];D=json.loads((ROOT/'inputs/sys.json').read_text());R=s.symbols(D['rvars']);a,b,c=s.symbols('a b c');mass=[a,a,b,b,c]
pairs=list(itertools.combinations(range(5),2));dist=dict(zip(pairs,R))
def r(i,j):return s.Integer(0) if i==j else dist[tuple(sorted((i,j)))]
def ac(i,j):
 return sum(mass[k]*(r(i,k)**-3-1)*(r(j,k)**2-r(i,k)**2-r(i,j)**2) for k in range(5) if k!=i)
def canonical(f):
 poly=s.Poly(s.expand(f),*R);terms=poly.terms();shift=tuple(min(m[j] for m,c in terms) for j in range(10))
 # Divide out distance monomial then normalize by the leading coefficient
 # in full Q[a,b,c,r] polynomial ring (a nonzero rational scalar only).
 g=s.expand(sum(coef*s.prod(x**(e-d) for x,e,d in zip(R,m,shift)) for m,coef in terms));p=s.Poly(g,*R,a,b,c,domain=s.QQ)
 return s.expand(g/p.LC())
expected=[];start=time.monotonic()
for i in range(5):
 for j in range(5):
  if i==j:continue
  den=s.prod(r(i,k)**3 for k in range(5) if k!=i)
  expected.append((f'AC({i+1},{j+1})',canonical(s.expand(ac(i,j)*den))))
for i,j in pairs:
 den=s.prod(r(i,k)**3 for k in range(5) if k!=i)*s.prod(r(j,k)**3 for k in range(5) if k!=j)
 expected.append((f'symmetric_AC({i+1},{j+1})',canonical(s.expand((ac(i,j)+ac(j,i))*den))))
for ids in itertools.combinations(range(5),4):
 mat=s.zeros(5)
 for k in range(1,5):mat[0,k]=mat[k,0]=1
 for i in range(4):
  for j in range(4):mat[i+1,j+1]=r(ids[i],ids[j])**2
 expected.append(('CM'+str(tuple(i+1 for i in ids)),canonical(mat.det(method='domain-ge'))))
available=list(expected);matches=[];failures=[]
for idx,f in enumerate(D['polys']):
 g=canonical(s.sympify(f));found=next((k for k,(label,h) in enumerate(available) if g==h),None)
 if found is None:failures.append(idx)
 else:matches.append({'stored_generator':idx,'formula':available.pop(found)[0]})
out={'matched':len(matches),'failures':failures,'unmatched_formulas':[label for label,g in available],'matches':matches,'elapsed_seconds':time.monotonic()-start,'normalization':'Distance monomial unit and nonzero rational scalar only; mass factors not removed.'}
(ROOT/'results/equation_verification.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='matches'});assert not failures and not available
