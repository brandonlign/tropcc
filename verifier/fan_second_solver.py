"""Replay fan completeness with Yices and standard SMT-LIB rational arithmetic.

Uses a Boolean cardinality encoding, independently of Z3's PB extension.
All six cases must be UNSAT; partial runs, unknown and timeouts fail.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,time
import verify as V
from fan_completeness import maximal_cones,check_symmetries,halfspaces

ROOT=Path(__file__).resolve().parents[1]

def integer(n):
    return str(n) if n>=0 else '(- %d)'%(-n)

def linear(row):
    terms=['(* %s w%d)'%(integer(v),i) for i,v in enumerate(row) if v]
    return '(+ %s)'%' '.join(terms) if len(terms)>1 else terms[0] if terms else '0'

def cardinality(terms):
    # At least one true; any true member requires a different true member.
    return ['(assert (or %s))'%' '.join(terms)]+[
        '(assert (or (not %s) %s))'%(term,' '.join(terms[:j]+terms[j+1:]))
        for j,term in enumerate(terms)]

def build_query(polynomials,inequalities):
    lines=['(set-logic QF_LRA)']+['(declare-fun w%d () Real)'%j for j in range(10)]
    for i,f in enumerate(polynomials):
        maximum='maximum%d'%i;lines.append('(declare-fun %s () Real)'%maximum);terms=[]
        for j,m in enumerate(f):
            weight=linear(m);name='top%d_%d'%(i,j);terms.append(name)
            lines += ['(assert (>= %s %s))'%(maximum,weight),
                      '(define-fun %s () Bool (= %s %s))'%(name,maximum,weight)]
        lines+=cardinality(terms)
    for cone in inequalities:
        conditions=['(%s %s 0)'%('=' if equal else '>=',linear(row)) for equal,row in cone]
        lines.append('(assert (not (and %s)))'%' '.join(conditions))
    lines += ['(assert (and (>= w%d (- 1)) (<= w%d 1)))'%(j,j) for j in range(10)]
    lines += ['(assert (or %s))'%' '.join('(= w%d %s)'%(j,integer(sign)) for j in range(10) for sign in (-1,1))]
    return '\n'.join(lines)+'\n'

def case_constraint(case):
    coordinate=(0,1,3)[case//2];sign=(-1,1)[case%2]
    lines=['(assert (= w%d %s))'%(coordinate,integer(sign))]
    if coordinate==0:
        lines += ['(assert (>= w1 w%d))'%j for j in (2,4,5)]
    elif coordinate==1:
        lines.append('(assert (>= w0 w7))')
    else:
        lines.append('(assert (>= w1 w2))')
    return '\n'.join(lines)+'\n(check-sat)\n'

def solve(query,timeout):
    start=time.monotonic()
    command=['yices-smt2','--timeout='+str(timeout)]
    run=subprocess.run(command,input=query,text=True,
                       capture_output=True,timeout=timeout+120,check=True)
    lines=run.stdout.strip().splitlines();assert len(lines)==1,run.stdout+run.stderr
    return lines[0],time.monotonic()-start

def controls():
    header='(set-logic QF_LRA)\n(declare-fun w0 () Real)\n(declare-fun w1 () Real)\n(declare-fun m () Real)\n'
    header+='(assert (and (>= m 0) (>= m w0) (>= m w1)))\n'
    header+='\n'.join(cardinality(['(= m 0)','(= m w0)','(= m w1)']))+'\n'
    rays=['(and (= w0 0) (<= w1 0))','(and (= w1 0) (<= w0 0))','(and (= w0 w1) (>= w0 0))']
    for n,expected in [(3,'unsat'),(2,'sat')]:
        result,_=solve(header+'(assert (not (or %s)))\n(check-sat)\n'%' '.join(rays[:n]),10)
        assert result==expected,(result,expected)
    print('PASS Yices complete/incomplete tropical-line controls',flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=int,choices=range(6))
    parser.add_argument('--timeout',type=int,default=3600)
    parser.add_argument('--halfspaces',type=Path,help='Development cache; full replay regenerates')
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args();controls()
    if args.self_test:
        return
    rays,cells=V.read_fan(ROOT/'data/fan.out');rv,polys=V.load_sys()
    cones=maximal_cones(ROOT/'data/fan.out');check_symmetries(rays,cells,rv,polys)
    inequalities=json.loads(args.halfspaces.read_text()) if args.halfspaces else halfspaces(rays,cones)
    assert len(inequalities)==len(cones)
    for cone,h in zip(cones,inequalities):
        for equal,row in h:
            assert all(V.dot(row,rays[r])==0 if equal else V.dot(row,rays[r])>=0 for r in cone)
    base=build_query(polys,inequalities);logs=ROOT/'logs';logs.mkdir(exist_ok=True)
    (logs/'fan_second_solver.smt2').write_text(base)
    print('Built standard SMT-LIB constraints',flush=True);outcomes=[]
    for case in range(6) if args.case is None else [args.case]:
        query=base+case_constraint(case);status,seconds=solve(query,args.timeout)
        record={'case':case,'coordinate':(0,1,3)[case//2],'sign':(-1,1)[case%2],
                'result':status,'seconds':seconds,'query_sha256':hashlib.sha256(query.encode()).hexdigest()}
        outcomes.append(record);print(record,flush=True)
    passed=args.case is None and all(r['result']=='unsat' for r in outcomes)
    report={'complete':passed,'solver':subprocess.check_output(['yices-smt2','--version'],text=True).strip(),
            'cases':outcomes,'fan_sha256':hashlib.sha256((ROOT/'data/fan.out').read_bytes()).hexdigest(),
            'system_sha256':hashlib.sha256((ROOT/'data/sys.json').read_bytes()).hexdigest(),
            'base_query_sha256':hashlib.sha256(base.encode()).hexdigest(),'halfspaces_rebuilt':args.halfspaces is None}
    name='yices_report'+('' if args.case is None else '_'+str(args.case))+'.json'
    (logs/name).write_text(json.dumps(report,indent=2)+'\n')
    if not passed:
        raise SystemExit(2)
    print('PASS independent Yices completeness: all six cases UNSAT',flush=True)

if __name__=='__main__':
    main()
