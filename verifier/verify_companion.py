"""Recheck the companion Laurent ledger without trusting .ok markers.

Reuse its independent Laurent identity checker, preserve its frozen reports,
then replay structural and transfer proofs and check the exact cell union.
"""
from pathlib import Path
import hashlib,json,os,subprocess,sys
import verify as V

PROJECT=Path(__file__).resolve().parents[1]
ROOT=PROJECT/'companion'
MASTER=ROOT/'checkers/verify_binomial.py'
def run_checker(script):
    # Its mathematical assertions are replayed, but runtime metadata must
    # not alter the frozen source/input manifest after verification.
    frozen={path:path.read_bytes() for path in (ROOT/'results').glob('*.json')}
    try:
        return subprocess.run([sys.executable,str(ROOT/'checkers'/script)],cwd=ROOT,
                              stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    finally:
        for path,data in frozen.items():
            if path.read_bytes()!=data:
                path.write_bytes(data)

names=['cell_binomial','equality','shifted','reduction','target16','target301',
       'pair','transfer','pair2','transfer2','target247','target279','target897','transfer897']
source=MASTER.read_text().split("out={'certificate_sha256'",1)[0]
union=set()
checked=[]
rays,cells=V.read_fan(PROJECT/'data/fan.out')
for name in names:
    sys.argv=[str(MASTER),'--input','results/'+name+'_certificates.json']
    env={'__file__':str(MASTER),'__name__':'__main__'}
    exec(compile(source,str(MASTER),'exec'),env)
    report=env['report']
    for certificate in report['certificates']:
        cell=cells[certificate['ray_index']]
        assert certificate['weight']==[sum(rays[r][j] for r in cell) for j in range(10)], 'certificate cell-weight mismatch'
    if name=='cell_binomial':
        sympy=env['s']; a,b,c=sympy.symbols('a b c')
        for i,f in enumerate(env['P']):
            monoms=sympy.Poly(f.as_expr(),a,b,c,*env['R']).monoms()
            assert all(sum(m[:3])==(1 if i<30 else 0) for m in monoms),'mass homogeneity'
        print('PASS mass homogeneity: c=1 reduction is valid',flush=True)
    union.update(c['ray_index'] for c in report['certificates'])
    checked.append({'batch':name,'count':len(report['certificates']),
                    'sha256':hashlib.sha256(env['source'].read_bytes()).hexdigest()})
    print('PASS Laurent batch',name,len(report['certificates']),flush=True)
for script in ['verify_structural.py','verify_ratio325_transfer.py',
               'audit_symmetry_closure.py','verify_equations.py','audit_cells.py']:
    print('RUN',script,flush=True)
    result=run_checker(script)
    print(result.stdout[-1500:],flush=True)
    assert result.returncode==0,script
union.update(json.loads((ROOT/'results/structural_verification.json').read_text())['positive_mass_cells'])
union.update(map(int,json.loads((ROOT/'results/ratio325_transfer_verification.json').read_text())['cells']))
closure=json.loads((ROOT/'results/symmetry_closure.json').read_text())
assert all(record['source_cell'] in union for record in closure['additional_generic_cells'].values()), 'unverified symmetry source'
union.update(map(int,json.loads((ROOT/'results/symmetry_closure.json').read_text())['additional_generic_cells']))
ids=set(json.loads((ROOT/'results/generic_ids.json').read_text())['generic'])
assert union==ids and len(ids)==10965,'ledger cell union mismatch'
for name in ('sys.json','fan.out'):
    assert (ROOT/'inputs'/name).read_bytes()==(PROJECT/'data'/name).read_bytes(),'input mismatch '+name
out={'verified':True,'generic_exclusions':len(ids),'batches':checked,
     'verifier_sha256':hashlib.sha256(MASTER.read_bytes()).hexdigest(),
     'inputs_match_main_project':True}
(PROJECT/'logs').mkdir(exist_ok=True)
(PROJECT/'logs/companion.json').write_text(json.dumps(out,indent=2)+'\n')
print('ALL PASS companion generic exclusions:',len(ids),flush=True)
