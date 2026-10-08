"""Verify paired-mass symmetries of the full input and locate cell-346 orbit."""
from pathlib import Path
import sympy as s,json,itertools,math,hashlib
ROOT=Path(__file__).resolve().parents[1];D=json.loads((ROOT/'inputs/sys.json').read_text());R=s.symbols(D['rvars']);a,b,c=s.symbols('a b c');pairs=list(itertools.combinations(range(5),2));polys=[s.sympify(f) for f in D['polys']]
def canonical(f):
 terms=s.Poly(s.expand(f),*R,a,b,c,domain=s.QQ);return s.expand(terms.as_expr()/terms.LC())
base={canonical(f) for f in polys};text=(ROOT/'inputs/fan.out').read_text()
def block(k):return text.split('\n'+k+'\n',1)[1].split('\n\n',1)[0].strip().splitlines()
rays=[tuple(int(l) for l in row.split('#')[0].split()) for row in block('RAYS')];cones=[tuple(map(int,row.split('#')[0].strip().strip('{}').split())) for row in block('CONES')];cones=[C for C in cones if C]
def primitive(w):
 g=math.gcd(*w);return tuple(i//g for i in w)
weights=[tuple(sum(rays[i][j] for i in C) for j in range(10)) for C in cones];lookup={}
for i,w in enumerate(weights):lookup.setdefault(primitive(w),[]).append(i)
w=weights[346];records=[]
for exchange in [False,True]:
 for flip1,flip2 in itertools.product([False,True],repeat=2):
  perm=[0,1,2,3,4]
  if flip1:perm[0],perm[1]=perm[1],perm[0]
  if flip2:perm[2],perm[3]=perm[3],perm[2]
  if exchange:perm=[{0:2,1:3,2:0,3:1,4:4}[i] for i in perm]
  mapping=[pairs.index(tuple(sorted((perm[i],perm[j])))) for i,j in pairs]
  sub={R[i]:R[j] for i,j in enumerate(mapping)}
  if exchange:sub.update({a:b,b:a})
  transformed={canonical(f.xreplace(sub)) for f in polys};assert transformed==base
  target=[0]*10
  for i,j in enumerate(mapping):target[j]=w[i]
  records.append({'body_permutation_0based':perm,'exchange_a_b':exchange,'full_polynomial_system_invariant':True,'target_weight':target,'matching_cells':lookup.get(primitive(target),[])})
ids=sorted({i for r in records for i in r['matching_cells']});out={'source_cell':346,'symmetries':records,'orbit_cells':ids,'input_sha256':{f:hashlib.sha256((ROOT/'inputs'/f).read_bytes()).hexdigest() for f in ['sys.json','fan.out']},'scope':'All positive masses preserved by parameter permutation; each listed cell has a positive-proportional transformed interior weight.'}
(ROOT/'results/symmetry_346.json').write_text(json.dumps(out,indent=2)+'\n');print('Verified symmetries:',len(records),'orbit cells:',ids)
