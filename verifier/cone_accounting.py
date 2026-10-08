"""Reconstruct disjoint proof-type counts from the verified coverage ledger.

This classifies evidence; it does not replace the full certificate replay.
"""
from pathlib import Path
from collections import Counter
import json
import verify as V

ROOT=Path(__file__).resolve().parents[1]

def main():
    _,cells=V.read_fan(ROOT/'data/fan.out')
    orbit={int(k):v for k,v in json.loads((ROOT/'verifier/orbits.json').read_text()).items()}
    sigma=json.loads((ROOT/'verifier/sigma.json').read_text())
    companion=set(json.loads((ROOT/'companion/results/generic_ids.json').read_text())['generic'])
    exact=ROOT/'data/exact'
    elimination={334,1767,49,2171,2174,2175,2345,218,2885,3494,3545}
    modular={81,316,121,50,206,220}
    def source(rep):
        for folder in ['cert','dcert','gcert','hcert','ecert']:
            path=exact/folder/f'{rep}.cof'
            if path.exists():
                with path.open() as stream:
                    header=[next(stream).strip() for _ in range(5)]
                kind='Laurent identities'
                if 'subset 35' in header:
                    kind='Scaling ray' if rep==79 else 'Single-term virial'
                return folder,kind
        for folder in ['dcert','gcert']:
            if (exact/folder/f'{rep}.res').exists():
                return folder,'Resultant specialization'
        if (exact/'split'/f'{rep}.split').exists():
            return 'hcert','Split trees'
        if rep in elimination:
            return 'elim','Elimination'
        if rep in modular:
            return 'lemS','Modular stars'
    directories=Counter(); kinds=Counter(); residual=[]; assignments={}; representatives={}
    for k in range(len(cells)):
        if k in companion:
            directories['companion']+=1;kinds['Laurent/structural companion']+=1
            assignments[str(k)]={'kind':'Laurent/structural companion','source_cell':k}
            continue
        rep=orbit[k]; found=source(rep); transport=False
        if found is None:
            rep=orbit[sigma[k]]; transported=source(rep)
            if transported and (exact/transported[0]/f'{rep}.cof').exists():
                found='sigma','Additional symmetry transport'
                transport=True
        if found:
            directories[found[0]]+=1;kinds[found[1]]+=1
            assignments[str(k)]={'kind':found[1],'representative':rep,'mass_exchange':transport}
            representatives.setdefault(found[1],{}).setdefault(str(rep),[]).append(k)
        else:
            residual.append(k)
    ledger=json.loads((ROOT/'verifier/coverage_report.json').read_text())
    assert dict(directories)==ledger['covered']
    assert residual==ledger['residual_cells']
    report={'counts':dict(kinds),'excluded':sum(kinds.values()),'residual':len(residual),
            'total':len(cells),'directory_counts':dict(directories),
            'representatives':representatives,'cell_assignments':assignments}
    (ROOT/'data/cone_accounting.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS disjoint accounting:',dict(kinds))
    for kind,reps in representatives.items():
        print(kind, len(reps),'representatives;',
              ', '.join(str(k)+':'+str(len(v)) for k,v in reps.items()) if len(reps)<30 else 'see JSON')

if __name__=='__main__':
    main()
