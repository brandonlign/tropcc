"""Verified finite group closure of existing direct/structural exclusions.
Positive-proportional interior weights transfer initial-system exclusion exactly.
"""
from pathlib import Path
import json,sympy as s,itertools,math,hashlib,sys
ROOT=Path(__file__).resolve().parents[1];D=json.loads((ROOT/'inputs/sys.json').read_text());R=s.symbols(D['rvars']);a,b,c=s.symbols('a b c');pairs=list(itertools.combinations(range(5),2));polys=[s.sympify(f) for f in D['polys']]
def canon(f):
 p=s.Poly(s.expand(f),*R,a,b,c,domain=s.QQ);return s.expand(p.as_expr()/p.LC())
base={canon(f) for f in polys};text=(ROOT/'inputs/fan.out').read_text()
def block(k):return text.split('\n'+k+'\n',1)[1].split('\n\n',1)[0].strip().splitlines()
rays=[tuple(map(int,row.split('#')[0].split())) for row in block('RAYS')];cones=[tuple(map(int,row.split('#')[0].strip().strip('{}').split())) for row in block('CONES')];cones=[C for C in cones if C];weights=[tuple(sum(rays[i][j] for i in C) for j in range(10)) for C in cones]
def primitive(w):return tuple(x//math.gcd(*w) for x in w)
lookup={}
for ci,w in enumerate(weights):lookup.setdefault(primitive(w),[]).append(ci)
manifest=json.loads((ROOT/'results/batches.json').read_text());generic={};uniform={};dependencies=['inputs/sys.json','inputs/fan.out','checkers/symmetry_closure.py','results/structural_verification.json','results/ratio325_transfer_verification.json']
batches=list(manifest['batches'])
for name in sys.argv[1:]:
 path=ROOT/'results'/f'{name}_certificates.json'
 if str(path.relative_to(ROOT)) not in {b['path'] for b in batches}:batches.append({'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
for batch in batches:
 path=ROOT/batch['path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==batch['sha256'];dependencies.append(batch['path'])
 d=json.loads(path.read_text())
 for cert in d['certificates']:
  ci=cert['ray_index'];generic[ci]={'kind':'direct','factor':cert['parameter_factor'],'batch':batch['path']}
  coeff=s.Poly(s.sympify(cert['parameter_factor']),a,b,c,domain=s.QQ).coeffs()
  if all(v>0 for v in coeff) or all(v<0 for v in coeff):uniform[ci]=generic[ci]
for ci in json.loads((ROOT/'results/structural_verification.json').read_text())['positive_mass_cells']:generic[ci]=uniform[ci]={'kind':'structural'}
for key,conditions in json.loads((ROOT/'results/ratio325_transfer_verification.json').read_text())['cells'].items():generic.setdefault(int(key),{'kind':'ratio325','necessary_positive_mass_equalities':conditions})
initial_generic=set(generic);initial_uniform=set(uniform);maps=[]
for exchange in [False,True]:
 for f1,f2 in itertools.product([False,True],repeat=2):
  perm=[0,1,2,3,4]
  if f1:perm[0],perm[1]=perm[1],perm[0]
  if f2:perm[2],perm[3]=perm[3],perm[2]
  if exchange:perm=[{0:2,1:3,2:0,3:1,4:4}[i] for i in perm]
  mapping=[pairs.index(tuple(sorted((perm[i],perm[j])))) for i,j in pairs];sub={R[i]:R[j] for i,j in enumerate(mapping)}
  if exchange:sub.update({a:b,b:a})
  assert {canon(f.xreplace(sub)) for f in polys}==base
  maps.append({'permutation':perm,'exchange_a_b':exchange,'distance_permutation':mapping})
new_generic={};new_uniform={};missing=0
for ci,proof in generic.items():
 for si,mapping in enumerate(maps):
  target=[0]*10
  for i,j in enumerate(mapping['distance_permutation']):target[j]=weights[ci][i]
  matches=lookup.get(primitive(target),[])
  if not matches:missing+=1
  for cj in matches:
   record={'source_cell':ci,'symmetry':si,'source_proof':proof}
   if cj not in initial_generic:new_generic.setdefault(cj,record)
   if ci in uniform and cj not in initial_uniform:new_uniform.setdefault(cj,record)
out={'symmetries':maps,'additional_generic_cells':new_generic,'additional_uniform_cells':new_uniform,'base_generic_count':len(initial_generic),'base_uniform_count':len(initial_uniform),'generic_union_count':len(initial_generic|set(new_generic)),'uniform_union_count':len(initial_uniform|set(new_uniform)),'unmatched_transformed_weights':missing,'dependency_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in dependencies},'scope':'Exact full-system symmetries and positively proportional interior weights. Generic conditions must be parameter-permuted, not discarded. Not a fan-completeness proof.'}
(ROOT/'results/symmetry_closure.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k not in ['symmetries','additional_generic_cells','additional_uniform_cells','dependency_sha256']});print('Additional generic',len(new_generic),'uniform',len(new_uniform))
