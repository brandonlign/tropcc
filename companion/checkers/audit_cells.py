"""Integer support audit: each listed cone lies in common max-normal cones.
Does not prove completeness of the list of cones.
"""
from pathlib import Path
import json, sympy as s
ROOT=Path(__file__).resolve().parents[1]
text=(ROOT/'inputs/fan.out').read_text()
def block(name): return text.split('\n'+name+'\n',1)[1].split('\n\n',1)[0].strip().splitlines()
rays=[tuple(map(int,l.split('#')[0].split())) for l in block('RAYS')]
cones=[tuple(map(int,l.split('#')[0].strip().strip('{}').split())) for l in block('CONES')]
cones=[c for c in cones if c]
D=json.loads((ROOT/'inputs/sys.json').read_text());R=s.symbols(D['rvars'])
mons=[s.Poly(s.sympify(f),*R).monoms() for f in D['polys']]
# Exact exposed supports on each ray; a nonempty intersection implies that
# the sum weight exposes precisely that intersection throughout relint(cone).
supports=[]; ray_not_tropical=[]
for ri,w in enumerate(rays):
    fs=[]
    for fi,ms in enumerate(mons):
        vals=[sum(a*b for a,b in zip(m,w)) for m in ms];top=max(vals)
        face=frozenset(i for i,v in enumerate(vals) if v==top);fs.append(face)
        if len(face)<2:ray_not_tropical.append([ri,fi])
    supports.append(fs)
failures=[]
for ci,c in enumerate(cones):
    for fi in range(len(mons)):
        common=set(supports[c[0]][fi])
        for ri in c[1:]:common.intersection_update(supports[ri][fi])
        if len(common)<2:failures.append([ci,fi,len(common)])
certs=json.loads((ROOT/'results/cell_binomial_certificates.json').read_text())['certificates']
weights_ok=all(tuple(z['weight'])==tuple(sum(rays[i][j] for i in cones[z['ray_index']]) for j in range(len(R))) for z in certs)
report={'ray_count':len(rays),'nonzero_cells':len(cones),'polynomials':len(mons),
        'ray_support_failures':ray_not_tropical,'cell_support_failures':failures,
        'certificate_weights_match_cell_sums':weights_ok,
        'lineality_dimension_recorded':block('LINEALITY_DIM'),
        'scope':'Exact max-normal support compatibility, not completeness or irredundancy of fan.'}
(ROOT/'results/cell_support_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:(len(v) if k.endswith('failures') else v) for k,v in report.items()},indent=2))
assert not failures and not ray_not_tropical and weights_ok
