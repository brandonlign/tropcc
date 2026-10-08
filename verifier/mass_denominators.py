"""Evaluate recorded cofactor denominators at a named rational mass pair.

This is a necessary certificate-specialization check, not certification of
all structural and modular exclusions at that point.
"""
from pathlib import Path
import argparse,json
import flint
import verify as V

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--a',type=int,default=3);parser.add_argument('--b',type=int,default=7)
    args=parser.parse_args();ctx=flint.fmpq_mpoly_ctx.get(('a','b'),'lex')
    values={};zero_files={};terms=0;files=0
    for directory in ['cert','dcert','gcert','hcert','ecert','split']:
        for path in sorted((ROOT/'data/exact'/directory).glob('*')):
            if path.suffix not in ('.cof','.split'):
                continue
            files+=1;zeros=set()
            with path.open() as stream:
                for line in stream:
                    if '|' not in line:
                        continue
                    parts=line.strip().split('|'); assert len(parts)==3
                    denominator=parts[1];terms+=1
                    if denominator not in values:
                        f=V.parse_poly(denominator,ctx,['a','b'])
                        assert all(c.denominator==1 for c in f.to_dict().values())
                        values[denominator]=int(f(args.a,args.b))
                    if values[denominator]==0:
                        zeros.add(denominator)
            if zeros:
                zero_files[str(path.relative_to(ROOT))]=sorted(zeros)
    report={'mass_pair':[args.a,args.b],'files_scanned':files,'cofactor_terms':terms,
            'distinct_denominators':len(values),'zero_denominators':sum(v==0 for v in values.values()),
            'files_with_zero_denominators':zero_files,
            'scope':'Explicit cofactor and split-tree files only; not an explicit instance of the theorem.'}
    (ROOT/'data/mass_denominator_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS denominator scan:',{k:v for k,v in report.items() if k!='files_with_zero_denominators'})
    if zero_files:
        print('Vanishing files:',list(zero_files))

if __name__=='__main__':
    main()
