#!/usr/bin/env python3
"""Resultant-specialisation certificates for 2-variable cell systems (S12). No Singular, exact (python-flint).
Certificate file (dcert/<k>.res or gcert/<k>.res), lines:
   cell k / [system F+G] / divmono / subset i1,i2,i3 / vars y,x / spec a0,b0
Checks (steps 3,4 identical to verify.py: CELL_LEMMA support constancy + torus slice onto the fixed coords):
  g1,g2,g3 = divmono initial forms in Z[a,b][x,y].  A = Res_x(g1,g2)/y^k1, B = Res_x(g1,g3)/y^k2 (y-powers removed:
  units of the Laurent ring).  A,B lie in <g> (Laurent, over Q[a,b]).  At the rational point (a0,b0): deg_y A, deg_y B are
  preserved and gcd(A(a0,b0), B(a0,b0)) = 1 in Q[y], hence D := Res_y(A,B) in Q[a,b] has D(a0,b0) = Res_y(A0,B0) != 0.
  So D is a NONZERO element of Q[a,b] lying in <g>: the ideal is the unit ideal over Q(a,b) (and at every (a,b) with D != 0).
"""
import sys,os
import flint
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import verify as VF
def check(fn):
    rv,F=VF.load_sys(); rays,cells=VF.read_fan(os.path.join(VF.N5,'fan.out'))
    L=[l.strip() for l in open(fn).read().split('\n') if l.strip()]
    k=int(L[0].split()[1]); L=L[1:]
    aug=L[0]=='system F+G'
    if aug: L=L[1:]; rv,F=VF.augment(rv,F)
    assert L[0]=='divmono'; L=L[1:]
    S=[int(x) for x in L[0].split()[1].split(',')]; assert len(S)==3
    y,x=L[1].split()[1].split(','); keep=[y,x]; assert y!=x and y in rv and x in rv
    a0,b0=[flint.fmpq(*map(int,(t.split('/')+['1'])[:2])) for t in L[2].split()[1].split(',')]
    fixj=[j for j,v in enumerate(rv) if v not in keep]
    V=[rays[i]+[0]*(len(rv)-len(rays[i])) for i in cells[k]]; w=[sum(z) for z in zip(*V)]
    forms=[]
    for i in S:
        sup=list(F[i]); Fi=set(sup)
        for v in V:
            M=max(VF.dot(v,m) for m in sup); Fi&={m for m in sup if VF.dot(v,m)==M}
        Mw=max(VF.dot(w,m) for m in sup); Fw={m for m in sup if VF.dot(w,m)==Mw}
        if not Fi or Fi!=Fw: return k,False,'support not constant on cell for f%d'%i
        forms.append(sorted(Fw))
    rows=[[p-q for p,q in zip(f[0],m)] for f in forms for m in f[1:]]
    H=VF.nullspace_Q(rows,len(rv))
    if fixj:
        P=flint.fmpq_mat(len(H),len(fixj),[h[j] for h in H for j in fixj]) if H else None
        if P is None or P.rank()!=len(fixj): return k,False,'torus slice invalid'
    names=['a','b',y,x]; ctx=flint.fmpz_mpoly_ctx.get(tuple(names),'lex'); kpos=[rv.index(y),rv.index(x)]
    G=[]
    for Fw,i in zip(forms,S):
        d={}
        for m in Fw:
            for (ea,eb),c in F[i][m].items():
                e=(ea,eb,m[kpos[0]],m[kpos[1]]); d[e]=d.get(e,0)+c
        mn=[min(e[2+t] for e in d) for t in range(2)]
        G.append(ctx.from_dict({(e[0],e[1],e[2]-mn[0],e[3]-mn[1]):c for e,c in d.items() if c}))
    Y=ctx.gens()[2]
    def strip(p):
        assert not p.is_zero(), 'zero resultant'
        while True:
            q,r=divmod(p,Y)
            if r!=0: return p
            p=q
    A=strip(G[0].resultant(G[1],x)); B=strip(G[0].resultant(G[2],x))
    assert A.degrees()[3]==0 and B.degrees()[3]==0
    def spec(p):
        c=[flint.fmpq(0)]*(p.degrees()[2]+1)
        for (ea,eb,ey,ex),v in p.to_dict().items(): c[ey]+=flint.fmpq(int(v))*a0**ea*b0**eb
        return flint.fmpq_poly(c)
    A0,B0=spec(A),spec(B)
    if A0.degree()!=A.degrees()[2] or B0.degree()!=B.degrees()[2]: return k,False,'leading coefficient vanishes at spec point'
    g=A0.gcd(B0)
    if g.degree()!=0: return k,False,'specialised gcd nontrivial'
    return k,True,'resultant cert OK [degs %d,%d; spec %s,%s]%s'%(A0.degree(),B0.degree(),a0,b0,' [F+G]' if aug else '')
if __name__=='__main__':
    bad=0
    for fn in sys.argv[1:]:
        try: k,ok,msg=check(fn)
        except Exception as e: k,ok,msg=fn,False,'malformed: %r'%e
        print(k,'PASS' if ok else 'FAIL',msg,flush=True); bad+=not ok
    print('failures',bad); sys.exit(1 if bad else 0)
