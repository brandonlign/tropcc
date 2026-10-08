#!/usr/bin/env python3
"""S24: transport the elimination certificates (verify_334.py, verify_49.py) along sigma=(13)(24), a<->b.
For each base cell r of a group, the header (system, divmono, subset, vars) is transported by transport.py to the cell sigma(r);
the forms of sigma(r) are recomputed by verify.py (support constancy + torus slice) and compared with those of r:
generator by generator they must agree up to a nonzero rational factor after a<->b (variables correspond positionally).
Then the polynomial system of sigma(r) is the system of r with a,b exchanged, and the elimination proof (generic (a,b)) applies verbatim.
usage: elim_sigma.py  -> prints the sigma-cells covered, exit 0 iff all checks pass"""
import sys,os,json,tempfile
from fractions import Fraction
import verify as V, transport
HERE=os.path.dirname(os.path.abspath(__file__))
SIG=json.load(open(os.path.join(HERE,'sigma.json')))
GROUPS={'verify_334.py':(['system F+G+H'],'subset 1,2,4,5,6,13','vars r13,r14,r23,r24,r34,u',[334,1767]),
        'verify_49.py':(['system F+G+H'],'subset 2,3,8,10,11,18','vars r14,r15,r34,r35,r45,u',[49,2171,2174,2175])}
def forms(lines):
    r=V.verify(None,lines=lines+[''],ret_gens=True); assert len(r)==4, r
    return r[3][:-1]
def same_upto_scalar_swap(G1,G2,swap=True):
    import flint
    for g,h in zip(G1,G2):
        dg={((e[1],e[0])+e[2:] if swap else e):Fraction(int(c.p),int(c.q)) if hasattr(c,'p') else Fraction(c) for e,c in ((e,flint.fmpq(c)) for e,c in g.to_dict().items())}
        dh={e:Fraction(int(flint.fmpq(c).p),int(flint.fmpq(c).q)) for e,c in h.to_dict().items()}
        if set(dg)!=set(dh): return False
        k0=next(iter(dg)); l=dh[k0]/dg[k0]
        if l==0 or any(dh[k]!=l*dg[k] for k in dg): return False
    return True
covered=[]; bad=0
for sc,(sysl,sub,vs,cells) in GROUPS.items():
    for r in cells:
        hdr=['cell %d'%r]+sysl+['divmono',sub,vs]
        G1=forms(hdr)
        t=SIG[r]
        if t==r: print(sc,r,'is sigma-invariant (nothing to transport)'); continue
        fin=os.path.join(tempfile.gettempdir(),'esig_in.cof'); fout=os.path.join(tempfile.gettempdir(),'esig_out.cof')
        open(fin,'w').write('\n'.join(hdr+['rhs prod^1','#1'])+'\n')
        transport.run(fin,t,fout)
        h2=[l for l in open(fout).read().split('\n') if l.strip() and not l.startswith('rhs') and not l.startswith('#')]
        G2=forms(h2)
        ok=same_upto_scalar_swap(G1,G2) or same_upto_scalar_swap([g for g in G1],G2,swap=False)
        print(sc,r,'-> sigma cell',t,'header',h2[1:],'PASS' if ok else 'FAIL',flush=True)
        if ok: covered.append(t)
        else: bad+=1
print('sigma-covered cells',covered,'failures',bad)
sys.exit(1 if bad else 0)
