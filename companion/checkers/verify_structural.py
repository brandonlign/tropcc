"""Rerun complete structural-proof dependency chain before manifest integration."""
from pathlib import Path
import subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
names=['eliminate_ratio_346','audit_ratio_346','sturm_mass_factor','exceptional_branch_346','audit_exceptional_346','symmetry_346']
for name in names:subprocess.run([sys.executable,str(ROOT/'checkers'/f'{name}.py')],check=True,timeout=60)
d=json.loads((ROOT/'results/symmetry_346.json').read_text())
artifacts=[ROOT/'checkers'/f'{n}.py' for n in names]+[ROOT/'inputs/sys.json',ROOT/'inputs/fan.out']
out={'verified':True,'positive_mass_cells':d['orbit_cells'],'proof_type':'Exact elimination, positive-mass Sturm classification, exceptional-branch contradiction, polynomial-system symmetries','dependency_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in artifacts},'scope':'No assertion of complete fan coverage or overall finiteness.'}
(ROOT/'results/structural_verification.json').write_text(json.dumps(out,indent=2)+'\n');print('STRUCTURAL PROOF VERIFIED',d['orbit_cells'])
