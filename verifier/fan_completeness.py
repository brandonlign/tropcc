"""Check support-fan completeness with cddlib/GMP and exact rational SMT.

Requires cddexec_gmp and z3-solver. Exit 0 requires all six symmetry cases
UNSAT; UNKNOWN (including timeout) exits 2 and is never a proof result.
"""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import argparse, hashlib, json, math, subprocess, time
import verify as V

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def maximal_cones(path):
    lines = path.read_text().split('MAXIMAL_CONES_OF_CLOSURE\n',1)[1].splitlines()
    result=[]
    for line in lines:
        if not line.strip():
            break
        text=line.split('#')[0].strip()
        assert text.startswith('{') and text.endswith('}')
        result.append(tuple(map(int,text[1:-1].split())))
    return result

def permute(v, permutation):
    out=[0]*len(v)
    for i,j in enumerate(permutation):
        out[j]=v[i]
    return tuple(out)

def check_symmetries(rays, cones, rv, polynomials):
    pairs=[tuple(map(int,r[1:])) for r in rv]
    index={p:i for i,p in enumerate(pairs)}
    supports=Counter(frozenset(f) for f in polynomials)
    ray_cones={frozenset(tuple(rays[i]) for i in c) for c in cones}
    for relabel in [(2,1,3,4,5),(1,2,4,3,5),(3,4,1,2,5)]:
        permutation=[index[tuple(sorted((relabel[i-1],relabel[j-1])))] for i,j in pairs]
        transformed=Counter(frozenset(permute(m,permutation) for m in f) for f in polynomials)
        assert transformed == supports
        assert {frozenset(permute(rays[i],permutation) for i in c) for c in cones} == ray_cones
    # Under this order-eight action the distance-coordinate orbits have
    # representatives r12, r13, r15 (indices 0,1,3).
    assert rv == ['r12','r13','r14','r15','r23','r24','r25','r34','r35','r45']

def halfspaces(rays, cones):
    result=[]
    for k,cone in enumerate(cones):
        rows=['1 '+' '.join(['0']*10)]
        rows += ['0 '+' '.join(map(str,rays[i])) for i in cone]
        src='V-representation\nbegin\n%d 11 rational\n%s\nend\n'%(len(rows),'\n'.join(rows))
        run=subprocess.run(['cddexec_gmp','--rep'],input=src,text=True,capture_output=True,check=True)
        lines=run.stdout.splitlines(); start=lines.index('begin')
        count,n,kind=lines[start+1].split(); assert n=='11' and kind=='rational'
        equal=set()
        for line in lines[:start]:
            if line.startswith('linearity'):
                data=line.split(); equal={int(x)-1 for x in data[2:]}
                assert len(equal)==int(data[1])
        h=[]
        for j,line in enumerate(lines[start+2:start+2+int(count)]):
            row=list(map(Fraction,line.split())); assert len(row)==11
            if row[0]:
                assert row[0]>0 and not any(row[1:]) and j not in equal
                continue
            denominator=math.lcm(*(x.denominator for x in row))
            normal=[int(x*denominator) for x in row[1:]]
            divisor=math.gcd(*normal); assert divisor>0
            normal=[x//divisor for x in normal]
            for i in cone:
                value=sum(a*b for a,b in zip(normal,rays[i]))
                assert value==0 if j in equal else value>=0
            h.append((j in equal,normal))
        result.append(h)
        if k%250==0:
            print('Converted',k,'of',len(cones),'cones',flush=True)
    return result

def build_solver(polynomials, halfspace_cones, z3, timeout, algorithm):
    w=z3.Reals('w0 w1 w2 w3 w4 w5 w6 w7 w8 w9')
    solver=z3.Solver(); solver.set(timeout=timeout); solver.set('arith.solver',algorithm)
    for i,f in enumerate(polynomials):
        maximum=z3.Real('max%d'%i)
        weights=[z3.Sum([v*x for v,x in zip(e,w) if v]) for e in f]
        solver.add(*(maximum>=x for x in weights))
        solver.add(z3.AtLeast(*[maximum==x for x in weights],2))
    for cone in halfspace_cones:
        inside=[]
        for equal,normal in cone:
            expr=z3.Sum([v*x for v,x in zip(normal,w) if v])
            inside.append(expr==0 if equal else expr>=0)
        solver.add(z3.Not(z3.And(*inside)))
    solver.add(*[z3.And(x>=-1,x<=1) for x in w])
    solver.add(z3.Or(*[z3.Or(x==1,x==-1) for x in w]))
    return solver,w

def self_test(z3):
    # The tropical hypersurface of 1+x+y comprises three rays. Removing
    # one ray must produce a rational witness outside the supplied fan.
    x,y,m=z3.Reals('x y m')
    tropical=z3.And(m>=0,m>=x,m>=y,z3.AtLeast(m==0,m==x,m==y,2))
    rays=[z3.And(x==0,y<=0),z3.And(y==0,x<=0),z3.And(x==y,x>=0)]
    complete=z3.Solver(); complete.add(tropical,z3.Not(z3.Or(*rays)))
    assert complete.check()==z3.unsat
    incomplete=z3.Solver(); incomplete.add(tropical,z3.Not(z3.Or(*rays[:2])))
    assert incomplete.check()==z3.sat
    print('PASS control: complete toy fan UNSAT, missing-ray fan SAT',flush=True)

def main():
    import z3
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--timeout',type=int,default=600,help='Seconds per symmetry case')
    parser.add_argument('--case',type=int,choices=range(6),help='Run one case only; cannot certify full completeness')
    parser.add_argument('--symmetry-break',action='store_true',help='Impose additional valid stabilizer inequalities')
    parser.add_argument('--arithmetic-solver',type=int,choices=(2,6),default=2)
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    self_test(z3)
    if args.self_test:
        return
    rays,cells=V.read_fan(ROOT/'data/fan.out')
    cones=maximal_cones(ROOT/'data/fan.out')
    stored={frozenset(c) for c in cells}
    assert all(frozenset(c) in stored for c in cones)
    rv,polynomials=V.load_sys()
    check_symmetries(rays,cones,rv,polynomials)
    inequalities=halfspaces(rays,cones)
    solver,w=build_solver(polynomials,inequalities,z3,args.timeout*1000,args.arithmetic_solver)
    logs=ROOT/'logs';logs.mkdir(exist_ok=True)
    (logs/'fan_completeness.smt2').write_text(solver.to_smt2())
    outcomes=[]
    for case in range(6) if args.case is None else [args.case]:
        coordinate=(0,1,3)[case//2]; sign=(-1,1)[case%2]
        solver.push(); solver.add(w[coordinate]==sign)
        if args.symmetry_break:
            if coordinate==0:
                solver.add(w[1]>=w[2],w[1]>=w[4],w[1]>=w[5])
            elif coordinate==1:
                solver.add(w[0]>=w[7])
            else:
                solver.add(w[1]>=w[2])
        start=time.monotonic(); status=solver.check()
        record={'case':case,'coordinate':coordinate,'sign':sign,'result':str(status),
                'seconds':time.monotonic()-start}
        if status==z3.sat:
            record['witness']=[str(solver.model().eval(x)) for x in w]
        elif status==z3.unknown:
            record['reason']=solver.reason_unknown()
        outcomes.append(record);print(record,flush=True);solver.pop()
    passed=args.case is None and all(r['result']=='unsat' for r in outcomes)
    report={'complete':passed,'z3_version':z3.get_version_string(),
            'fan_sha256':digest(ROOT/'data/fan.out'),'system_sha256':digest(ROOT/'data/sys.json'),
            'maximal_cones':len(cones),'symmetry_break':args.symmetry_break,'cases':outcomes}
    (logs/'fan_completeness_report.json').write_text(json.dumps(report,indent=2)+'\n')
    if passed:
        print('PASS independent completeness: no rational weight outside the stored fan')
    else:
        raise SystemExit(2)

if __name__=='__main__':
    main()
