#!/usr/bin/env python3
"""Reconstruct sliced modular cell forms for B0 certificate generation.

This is a search model, not a complete cell-exclusion proof. The independent
complete replay is verify_lemS_rep.py, called by verify_lemS.py.
"""
import sys,os,flint
import verify as V, verify_lemS as LS
CELL=int(sys.argv[1]) if len(sys.argv)>1 and sys.argv[1].isdigit() else 81
P0,A0,B0=32003,3,7
OK=[]
def check(name,cond):
    print(('PASS ' if cond else 'FAIL ')+name,flush=True); OK.append(bool(cond))
    if not cond: sys.exit(1)
rv=LS.rv
VV=('x','y','S','T','W','s23','s24','s25','r12','lam','u')
C=flint.nmod_mpoly_ctx.get(VV,P0,'lex'); x,y,S,T,W,s23,s24,s25,r12,lam,u=C.gens()
a,b=C.constant(A0),C.constant(B0)
CFG={81:dict(slice='r13',MAP={'r14':'x','r15':'y','r34':'S','r35':'T','r45':'W','r23':'s23','r24':'s24','r25':'s25','r12':'r12'},
           lin=4,g0=0,sq={'S':34,'T':33,'W':32},DST=9,DSW=13,leaf=(27,28)),
     121:dict(slice='r12',MAP={'r24':'x','r25':'y','r14':'S','r15':'T','r45':'W','r13':'s23','r34':'s24','r35':'s25','r23':'r12'},
           lin=9,g0=5,sq={'S':34,'T':33,'W':30},DST=1,DSW=14,leaf=(22,23))}
CFG[206]=dict(slice='r12',MAP={'r24':'x','r25':'y','r14':'S','r15':'T','r45':'W','r13':'s23','r23':'s24','r34':'s25','r35':'r12'},lin=35,g0=5,sq={'S':34,'T':33,'W':30},DST=1,DSW=14,leaf=(22,25))
K=CFG[CELL]; MAP=K['MAP']
def load():
    fs=LS.forms(CELL); check('0a all 37 forms of F+G+H constant and U-free on cell %d'%CELL,len(fs)==37)
    # torus slice r13=1: homogeneity space of the forms must project onto the r13 coordinate (verify.py step 4)
    rows=[[p_-q_ for p_,q_ in zip(min(f),m)] for f in fs.values() for m in f]
    H=V.nullspace_Q(rows,len(rv)); j13=rv.index(K['slice'])
    check('0b torus slice %s=1 valid'%K['slice']+' (homogeneity weight with nonzero r13 entry)',any(h[j13]!=0 for h in H))
    G={}
    for i,f in fs.items():
        mins=[min(m[j] for m in f) for j in range(len(rv))]; d={}
        for m,cf in f.items():
            e=[0]*len(VV)
            for j,n in enumerate(rv):
                k=m[j]-mins[j]
                if not k or n in (K['slice'],'U'): continue
                if MAP[n] in ('s23','s24','s25'):
                    if k%2: return None
                    k//=2
                e[VV.index(MAP[n])]=k
            d[tuple(e)]=(d.get(tuple(e),0)+LS.spec(cf,P0,A0,B0))%P0
        G[i]=C.from_dict({e:c for e,c in d.items() if c})
    return G
G=load(); check('0c r23,r24,r25 occur only with even exponents (s=r^2, surjective on torus, p odd)',G is not None)
def sub(f,**kw):
    gs=list(C.gens())
    for k,v in kw.items(): gs[VV.index(k)]=v
    return f.compose(*gs)
def coef(f,v,k):
    j=VV.index(v); return C.from_dict({tuple(0 if i==j else z for i,z in enumerate(e)):c for e,c in f.to_dict().items() if e[j]==k})
def deg(f,v): return f.degrees()[VV.index(v)] if not f.is_zero() else -1
def prop(f,g):  # f = c*g for a nonzero constant c ?
    if g.is_zero() or f.is_zero(): return False
    c=f.leading_coefficient()*pow(int(g.leading_coefficient()),-1,P0)
    return f==g*c
X,Y=x**3,y**3
gl=G[K['lin']]; l23,l24,l25=[coef(gl,v,1) for v in ('s23','s24','s25')]
check('1 linear form: l23 s23 + l24 s24 + l25 s25 with constant l25 != 0',gl==l23*s23+l24*s24+l25*s25 and l25.is_constant() and not l25.is_zero() and l23.is_constant() and l24.is_constant())
il25=pow(int(l25.leading_coefficient()),-1,P0)
G1={i:sub(g,s25=-(l23*s23+l24*s24)*il25) for i,g in G.items()}
g0=G1[K['g0']]; c23,c24=coef(g0,'s23',1),coef(g0,'s24',1)
check('2a g0|_{s25} = c23 s23 + c24 s24 (c in x,y)',g0==c23*s23+c24*s24 and all(deg(c,v)<=0 for c in (c23,c24) for v in ('S','T','W','r12')))
P,Q=c24,-c23
S25p=-(l23*P+l24*Q)*il25   # s25 = lam*S25p
check('2b P,Q,s25-coefficient nonzero polynomials (main branch: s23,s24,s25 != 0 => P,Q,S25p != 0)',not P.is_zero() and not Q.is_zero() and not S25p.is_zero())
Zx=x*y*P*Q*S25p
B0eqs=[c23,c24]   # branch B0: c23=c24=0

def mnorm(f):
    if f.is_zero(): return f
    d=f.to_dict(); m=[min(e[i] for e in d) for i in range(len(VV))]
    g=C.from_dict({tuple(z-m[i] for i,z in enumerate(e)):c for e,c in d.items()})
    return g*pow(int(g.leading_coefficient()),-1,P0)
LF={}
for i,g in sorted(G1.items()):
    h=sub(g,s23=lam*P,s24=lam*Q)
    if h.is_zero() or deg(h,'r12')>0: continue
    ks=[k for k in range(deg(h,'lam')+1) if not coef(h,'lam',k).is_zero()]
    if len(ks)==1: LF[i]=mnorm(coef(h,'lam',ks[0]))
