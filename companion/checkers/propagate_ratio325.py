"""Transfer structural proof when its used normalized initials are identical.
Only matching raw supports of the same generators is used; no heuristic lifting.
"""
from pathlib import Path
import json,sympy as s,itertools,hashlib
ROOT=Path(__file__).resolve().parents[1];D=json.loads((ROOT/'inputs/sys.json').read_text());R=s.symbols(D['rvars']);a,b,c=s.symbols('a b c');polys=[s.sympify(f) for f in D['polys']]
def canon(f):
 p=s.Poly(f,*R,a,b,c,domain=s.QQ);return s.expand(p.as_expr()/p.LC())
lookup={canon(f):i for i,f in enumerate(polys)};terms=[list(s.Poly(f,*R).monoms()) for f in polys]
text=(ROOT/'inputs/fan.out').read_text()
def block(k):return text.split('\n'+k+'\n',1)[1].split('\n\n',1)[0].strip().splitlines()
rays=[tuple(map(int,row.split('#')[0].split())) for row in block('RAYS')];cones=[tuple(map(int,row.split('#')[0].strip().strip('{}').split())) for row in block('CONES')];cones=[C for C in cones if C]
weights=[tuple(sum(rays[i][j] for i in C) for j in range(10)) for C in cones]
def support(i,w):
 scores=[sum(x*y for x,y in zip(m,w)) for m in terms[i]];top=max(scores);return tuple(k for k,v in enumerate(scores) if v==top)
used=[1,2,8];pairs=list(itertools.combinations(range(5),2));templates=[]
for row in json.loads((ROOT/'results/symmetry_346.json').read_text())['symmetries']:
 perm=row['body_permutation_0based'];sub={R[i]:R[pairs.index(tuple(sorted((perm[j],perm[k]))))] for i,(j,k) in enumerate(pairs)}
 if row['exchange_a_b']:sub.update({a:b,b:a})
 assert {canon(f.xreplace(sub)) for f in polys}==set(lookup)
 indices=[lookup[canon(polys[i].xreplace(sub))] for i in used];w=[0]*10
 for i,(j,k) in enumerate(pairs):w[pairs.index(tuple(sorted((perm[j],perm[k]))))]=weights[325][i]
 templates.append([(i,support(i,w)) for i in indices])
matched={}
for ci,w in enumerate(weights):
 cache={}
 for ti,template in enumerate(templates):
  ok=True
  for i,expected in template:
   if i not in cache:cache[i]=support(i,w)
   if cache[i]!=expected:ok=False;break
  if ok:matched.setdefault(ci,[]).append(ti)
out={'source_cell':325,'used_original_generators':used,'templates':templates,'matching_cells':matched,'count':len(matched),'method':'Identical raw initial supports for every used equation under a verified polynomial symmetry; original coefficients unchanged.','input_sha256':{f:hashlib.sha256((ROOT/'inputs'/f).read_bytes()).hexdigest() for f in ['sys.json','fan.out']}}
out['template_mass_conditions']=['a=4*b' if r['exchange_a_b'] else 'b=4*a' for r in json.loads((ROOT/'results/symmetry_346.json').read_text())['symmetries']]
(ROOT/'results/ratio325_propagation.json').write_text(json.dumps(out,indent=2)+'\n');print('Matched cells:',len(matched))
