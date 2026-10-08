"""Independent coefficient-level audit of every reported conditional transfer."""
from pathlib import Path
import json,sympy as s,runpy,itertools,hashlib
ROOT=Path(__file__).resolve().parents[1];ctx=runpy.run_path(str(ROOT/'checkers/ratio_325.py'));R=ctx['R'];a,b,c=ctx['a'],ctx['b'],ctx['c'];D=ctx['D'];rays=ctx['rays'];cones=ctx['cones'];initial=ctx['initial'];pairs=list(itertools.combinations(range(5),2));report=json.loads((ROOT/'results/ratio325_propagation.json').read_text());syms=json.loads((ROOT/'results/symmetry_346.json').read_text())['symmetries'];polys=[s.sympify(f) for f in D['polys']]
def canonical(f):
 terms=s.Poly(s.expand(f),*R).terms();shift=tuple(min(m[j] for m,coef in terms) for j in range(10));f=s.expand(sum(coef*s.prod(x**(e-d) for x,e,d in zip(R,m,shift)) for m,coef in terms));P=s.Poly(f,*R,a,b,c,domain=s.QQ);return s.expand(f/P.LC())
expected=[]
for row in syms:
 perm=row['body_permutation_0based'];sub={R[i]:R[pairs.index(tuple(sorted((perm[j],perm[k]))))] for i,(j,k) in enumerate(pairs)}
 if row['exchange_a_b']:sub.update({a:b,b:a})
 expected.append([canonical(initial(i).xreplace(sub)) for i in [1,2,8]])
cache={};checked=0;conditions={}
for ci,tis in report['matching_cells'].items():
 w=tuple(sum(rays[i][j] for i in cones[int(ci)]) for j in range(10));cs=set()
 for ti in tis:
  for position,(fi,unused_support) in enumerate(report['templates'][ti]):
   terms=s.Poly(polys[fi],*R).terms();scores=[sum(x*y for x,y in zip(m,w)) for m,coef in terms];top=max(scores);terms=[term for term,score in zip(terms,scores) if score==top];key=(fi,tuple(m for m,coef in terms))
   if key not in cache:cache[key]=canonical(sum(coef*s.prod(x**e for x,e in zip(R,m)) for m,coef in terms))
   assert cache[key]==expected[ti][position]
  cs.add('a=4*b' if syms[ti]['exchange_a_b'] else 'b=4*a');checked+=1
 conditions[ci]=sorted(cs)
paths=['inputs/sys.json','inputs/fan.out','checkers/ratio_325.py','checkers/propagate_ratio325.py','checkers/verify_ratio325_transfer.py','results/ratio325_propagation.json','results/symmetry_346.json']
out={'verified':True,'cells':conditions,'verified_template_matches':checked,'dependency_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},'scope':'Each cell excluded for positive masses away from its stated necessary mass equality; generic complex exclusion follows from the source mass polynomial. Not uniform exclusion on the equality.'}
(ROOT/'results/ratio325_transfer_verification.json').write_text(json.dumps(out,indent=2)+'\n');print('Verified coefficient-level matches:',checked,'cells:',len(conditions))
