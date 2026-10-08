"""Independent weight and proof-source audit of newly transferred exclusions."""
from pathlib import Path
import json,hashlib,math
ROOT=Path(__file__).resolve().parents[1];report=json.loads((ROOT/'results/symmetry_closure.json').read_text())
for filename,digest in report['dependency_sha256'].items():assert hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()==digest
text=(ROOT/'inputs/fan.out').read_text()
def block(k):return text.split('\n'+k+'\n',1)[1].split('\n\n',1)[0].strip().splitlines()
rays=[tuple(map(int,row.split('#')[0].split())) for row in block('RAYS')];cones=[tuple(map(int,row.split('#')[0].strip().strip('{}').split())) for row in block('CONES')];cones=[C for C in cones if C]
def weight(ci):return tuple(sum(rays[i][j] for i in cones[ci]) for j in range(10))
checks=0
for category in ['additional_generic_cells','additional_uniform_cells']:
 for target,record in report[category].items():
  source=record['source_cell'];mapping=report['symmetries'][record['symmetry']]['distance_permutation'];w=weight(source);v=weight(int(target));g1=math.gcd(*w);g2=math.gcd(*v)
  assert all(w[i]*g2==v[j]*g1 for i,j in enumerate(mapping));assert g1>0 and g2>0
  # Every source is a base proof, not a circular symmetry transfer.
  assert record['source_proof']['kind'] in ['direct','structural','ratio325']
  if category=='additional_uniform_cells':assert record['source_proof']['kind'] in ['direct','structural']
  checks+=1
out={'verified_records':checks,'closure_sha256':hashlib.sha256((ROOT/'results/symmetry_closure.json').read_bytes()).hexdigest(),'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Checks content provenance, positive proportionality and noncircular source types; polynomial-system symmetry itself is verified by symmetry_closure.py.'}
(ROOT/'results/symmetry_closure_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
