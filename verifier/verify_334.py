#!/usr/bin/env python3
"""Elimination certificate for the 334-group (cells 334 and 1767; S19-S24). No Singular; exact arithmetic (python-flint).
Claim: for generic masses the 6 initial forms g1..g6 (subset 1,2,4,5,6,13 of cell 334 [F+G+H] and of cell 1767 [F],
divmono, keep vars r13,r14,r23,r24,r34) have NO common zero on the complex torus. Equivalently the ideal is the unit ideal
over Q(a,b)(Laurent).  The proof is pointwise: K = algebraic closure of Q(a) (a transcendental, b=1 by the mass torus),
z in (K^*)^5 a hypothetical common zero; every derived polynomial below vanishes at z (justification in comments); at the
end every branch reaches a contradiction.  Steps:
  0. forms via verify.py (CELL_LEMMA support constancy + torus slice), identical for 334 and 1767; mass-torus weight with w_b!=0 => b=1.
  1. g3,g4,g5 are linear in X=r13^2,Y=r14^2,S=r34^2; Cramer: D*r13^2-Xn etc. are explicit combinations of g3,g4,g5; D is a monomial.
  2. replay of the division-free elimination (S20) -> E1..E14 in Z[a][r23,r24], all vanishing at z.
  3. branch r23=r24 (diagonal, S22): r14=+-r13 from g4-g5; resultant chain -> nonzero R(a) in Z[a].
  4. branch r23!=r24 (S24): symmetry (r13 r14)(r23 r24); symmetric/antisymmetric parts in s=r23+r24, p=r23*r24;
     Cramer-Delta factors F1..F4 split off; leaf and F-branches by gcds of resultants over Z[a,p]; the only common zero
     p=3/(a^2-1), s=-3a/(a^2-1) has Xq*Yq*Sq=0, i.e. r13 r14 r34=0: not on the torus.
Usage: verify_334.py   (prints PASS/FAIL per step, exit code 0 iff all PASS)
"""
import sys,os,time
import flint
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import verify as VF
OK=[]
def check(name,cond):
    print(('PASS ' if cond else 'FAIL ')+name,flush=True); OK.append(bool(cond))
    if not cond: print('failures 1'); sys.exit(1)

# ---------- 0. the forms ----------
V=('a','b','r13','r14','r34','r23','r24'); ctx=flint.fmpz_mpoly_ctx.get(V,'lex')
a,b,r13,r14,r34,r23,r24=ctx.gens(); zero=0*a; one=zero+1
def getforms(k,sysl):
    L=['cell %d'%k]+sysl+['divmono','subset 1,2,4,5,6,13','vars r13,r14,r23,r24,r34,u','']
    r=VF.verify(None,lines=L,ret_gens=True)
    assert len(r)==4, r
    _,c0,names,G=r
    assert names==['a','b','r13','r14','r23','r24','r34','u']
    out=[]
    for g in G[:6]:
        d={}
        for e,v in g.to_dict().items():
            v=flint.fmpq(v); assert v.q==1 and e[7]==0
            d[(e[0],e[1],e[2],e[3],e[6],e[4],e[5])]=int(v.p)
        out.append(d)
    return out
F334=getforms(334,['system F+G+H']); F1767=getforms(1767,[])
check('0a forms of 334 (F+G+H) and 1767 (F) on subset 1,2,4,5,6,13 coincide',F334==F1767)
# mass torus: weights (w_a,w_b,w_r13,..) making every g_i quasi-homogeneous; need one with w_b != 0
rows=[[x-y for x,y in zip(es[0],e)] for d in F334 for es in [sorted(d)] for e in es[1:]]
Hn=VF.nullspace_Q(rows,7)
check('0b mass torus: a quasi-homogeneity weight with w_b != 0 exists (b=1 WLOG)',any(h[1]!=0 for h in Hn))
def b1(d):  # set b=1
    o={}
    for e,c in d.items(): k=(e[0],0)+e[2:]; o[k]=o.get(k,0)+c
    return ctx.from_dict({k:c for k,c in o.items() if c})
G=[b1(d) for d in F334]
if os.environ.get('CORRUPT'):   # corruption test: perturb one form, the verifier must FAIL
    j=int(os.environ['CORRUPT'])
    if j<6: G[j]=G[j]+a*r24**3*r34**(j==5)
    elif j<9: G[5]=G[5]+(j-5)*r14**3*r34**2*r24   # change one integer coefficient of g6 (still overdetermined: may PASS)
    else: G[5]=r34*G[0]+r24*G[3]   # g6 in <g1..g5>: g1..g5 has torus zeros, must FAIL

# ---------- 1. Cramer for the squares ----------
def split(g,vs):
    d={}
    for e,c in g.to_dict().items():
        k=tuple(e[i] for i in vs); r=list(e)
        for i in vs: r[i]=0
        d.setdefault(k,{}); d[k][tuple(r)]=c
    return {k:ctx.from_dict(v) for k,v in d.items()}
keys=[(2,0,0),(0,2,0),(0,0,2)]
Ls=[split(G[i],(2,3,4)) for i in (2,3,4)]
check('1a g3,g4,g5 affine-linear in r13^2,r14^2,r34^2',all(set(l)<=set(keys+[(0,0,0)]) for l in Ls))
M=[[l.get(k,zero) for k in keys] for l in Ls]; cv=[-l.get((0,0,0),zero) for l in Ls]
def det3(M): return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
D=det3(M)
num=[det3([[cv[i] if jj==j else M[i][jj] for jj in range(3)] for i in range(3)]) for j in range(3)]
Xn,Yn,Sn=num
def minor(i,j): return [[M[r][c] for c in range(3) if c!=j] for r in range(3) if r!=i]
adj=[[(-1)**(i+j)*(lambda m:m[0][0]*m[1][1]-m[0][1]*m[1][0])(minor(j,i)) for j in range(3)] for i in range(3)]
sq=[r13**2,r14**2,r34**2]
check('1b D*r^2 - num = sum adj*g (Cramer identities)',all(D*sq[j]-num[j]==sum((adj[j][i]*G[2+i] for i in range(3)),zero) for j in range(3)))
check('1c D is a nonzero monomial (unit on the torus)',len(D)==1)

# ---------- 2. division-free elimination (replay of S20, exact) ----------
# At z: D*r13^2=Xn, D*r14^2=Yn, D*r34^2=Sn (step 1).  red(g) = D^K g with all squares replaced => red(g)(z) = D^K g(z).
# strip() divides by the monomial*integer content: a and the r's are nonzero at z.  So every poly below vanishes at z.
def red(g):
    terms=[]
    for e,c in g.to_dict().items():
        h=(e[2]//2,e[3]//2,e[4]//2); m=ctx.from_dict({(e[0],e[1],e[2]%2,e[3]%2,e[4]%2,e[5],e[6]):c}); terms.append((h,m))
    K=max(sum(h) for h,_ in terms)
    return sum((m*Xn**h[0]*Yn**h[1]*Sn**h[2]*D**(K-sum(h)) for h,m in terms),zero)
def strip(p): return p//p.term_content() if not p.is_zero() else p
def coef(p,var,k):
    return ctx.from_dict({tuple(0 if i==var else x for i,x in enumerate(e)):c for e,c in p.to_dict().items() if e[var]==k})
t0=time.time()
R1,R2,R6=[strip(red(G[i])) for i in (0,1,5)]
A1,B1,C1=coef(coef(R1,2,1),3,1),coef(coef(R1,2,1),3,0),coef(coef(R1,2,0),3,1)
A2,C2,B2=coef(coef(R2,3,1),2,1),coef(coef(R2,3,1),2,0),coef(coef(R2,3,0),2,1)
A6,B6,C6=coef(coef(R6,4,1),3,1),coef(coef(R6,4,1),3,0),coef(coef(R6,4,0),3,1)
check('2a R1=r13(A1 r14+B1)+C1 r14, R2=r14(A2 r13+C2)+B2 r13, R6=r34(A6 r14+B6)+C6 r14 (coefficients free of r13,r14,r34)',
      R1==r13*(A1*r14+B1)+C1*r14 and R2==r14*(A2*r13+C2)+B2*r13 and R6==r34*(A6*r14+B6)+C6*r14
      and all(all(e[2]==e[3]==e[4]==0 for e in q.to_dict()) for q in (A1,B1,C1,A2,B2,C2,A6,B6,C6)))
# R1(z)=0 => r13(A1r14+B1) = -C1 r14 => C1^2 r14^2 - r13^2(A1 r14+B1)^2 = 0 at z (same for R2, R6)
L1=strip(red(C1**2*r14**2-r13**2*(A1*r14+B1)**2)); L2=strip(red(B2**2*r13**2-r14**2*(A2*r13+C2)**2))
L6=strip(red(C6**2*r14**2-r34**2*(A6*r14+B6)**2))
al,be=coef(L1,3,1),coef(L1,3,0); ga,de=coef(L2,2,1),coef(L2,2,0); al6,be6=coef(L6,3,1),coef(L6,3,0)
check('2b L1=al r14+be, L2=ga r13+de, L6=al6 r14+be6 with coefficients in Z[a][r23,r24]',
      L1==al*r14+be and L2==ga*r13+de and L6==al6*r14+be6
      and all(all(e[2]==e[3]==e[4]==0 for e in q.to_dict()) for q in (al,be,ga,de,al6,be6)))
# at z: al r14=-be, ga r13=-de, al6 r14=-be6 and D r14^2=Yn, D r13^2=Xn.  Hence (squaring / eliminating r13,r14):
E={}
E[1]=strip(be**2*D-al**2*Yn); E[2]=strip(de**2*D-ga**2*Xn); E[3]=strip(al6*be-al*be6)
E[4]=strip(A1*be*de-B1*al*de-C1*ga*be)          # al*ga*R1 at z
E[5]=strip(A2*de*be-C2*be*ga-B2*de*al)          # al*ga*R2 at z
# R1: r14(A1 r13+C1) = -B1 r13 ; times ga: r14(C1 ga - A1 de) = B1 de... squared with D r14^2 = Yn (and similarly):
E[6]=strip((B1*de)**2*D-(C1*ga-A1*de)**2*Yn); E[7]=strip((B2*de)**2*D-(C2*ga-A2*de)**2*Yn)
E[8]=strip((C1*be)**2*D-(B1*al-A1*be)**2*Xn); E[9]=strip((C2*be)**2*D-(B2*al-A2*be)**2*Xn)
def strippow(P,v):
    while True:
        q,r=divmod(P,v)
        if not r.is_zero(): return P
        P=q
# Res_{r13}(R1,R2) is in the ideal (R1,R2) => vanishes at z; then reduce squares and r-powers as before
M1=strip(red(strippow(R1.resultant(R2,'r13'),r14))); M2=strip(red(strippow(R1.resultant(R2,'r14'),r13)))
mu,nu=coef(M1,3,1),coef(M1,3,0); m2,n2=coef(M2,2,1),coef(M2,2,0)
check('2c M1=mu r14+nu, M2=m2 r13+n2 (coefficients in Z[a][r23,r24])',M1==mu*r14+nu and M2==m2*r13+n2
      and all(all(e[2]==e[3]==e[4]==0 for e in q.to_dict()) for q in (mu,nu,m2,n2)))
E[10]=strip(nu**2*D-mu**2*Yn); E[11]=strip(mu*be-al*nu); E[12]=strip(mu*be6-al6*nu); E[13]=strip(n2**2*D-m2**2*Xn); E[14]=strip(m2*de-ga*n2)
check('2d E1..E14 nonzero, in Z[a][r23,r24] (%.0fs)'%(time.time()-t0),
      all(not e.is_zero() and all(x[2]==x[3]==x[4]==0 for x in e.to_dict()) for e in E.values()))

# ---------- 3. diagonal branch r23=r24 ----------
cd=flint.fmpz_mpoly_ctx.get(('a','x','t','s'),'lex'); ax,x,t,s=cd.gens()
def diag(g,eps):  # r13=x, r14=eps*x, r23=r24=t, r34=s
    o=0*ax
    for e,c in g.to_dict().items(): o+=c*ax**e[0]*x**e[2]*(eps*x)**e[3]*s**e[4]*t**e[5]*t**e[6]
    return o
def stripd(p):
    tc=p.term_content(); return p//tc
dd=G[3]-G[4]
dsub=sum((c*a**e[0]*r13**e[2]*r14**e[3]*r34**e[4]*r23**(e[5]+e[6]) for e,c in dd.to_dict().items()),zero)   # r24 -> r23
check('3a (g4-g5)|_{r24=r23} = +-a r23^4 (r13-r14)(r13+r14)  => r14 = +-r13',dsub in (a*r23**4*(r13**2-r14**2),-a*r23**4*(r13**2-r14**2)))
Hm=[stripd(diag(g,-1)) for g in G]
check('3b eps=-1: g1/x - g2/x = +-2 r34^2 => r34=0, impossible',(Hm[0]-Hm[1]) in (2*s**2,-2*s**2) or (Hm[0]+Hm[1]) in (2*s**2,-2*s**2))
Hp=[stripd(diag(g,1)) for g in G]
def norm_s2(p,sg):   # normalise sign so that coefficient of s^2 (s-free part) is sg
    c=p.coefficient(p.monoms().index((0,0,0,2))) if (0,0,0,2) in p.monoms() else None
    return p if c==sg else -p
A_=norm_s2(Hp[0],-1); B_=Hp[2]; C_=norm_s2(Hp[3],1); D_=Hp[5]
check('3c A+C free of s (A=g1, C=g4 after r-monomial removal)',all(e[3]==0 for e in (A_+C_).monoms()))
t0=time.time()
Q1=stripd(D_.resultant(C_,'s')); Q2=stripd((A_+C_).resultant(B_,'x')); Q3=stripd(Q1.resultant(B_,'x'))
Rd=Q2.resultant(Q3,'t')
check('3d resultant chain R=Res_t(Res_x(A+C,B),Res_x(Res_s(g6,C),B)) is a nonzero element of Z[a] (%.0fs)'%(time.time()-t0),
      (not Rd.is_zero()) and all(e[1]==e[2]==e[3]==0 for e in Rd.monoms()))

# ---------- 4. off-diagonal branch r23 != r24 ----------
def sw(g):  # (r13 r14)(r23 r24)
    return sum((c*a**e[0]*r13**e[3]*r14**e[2]*r34**e[4]*r23**e[6]*r24**e[5] for e,c in g.to_dict().items()),zero)
SG=[sw(g) for g in G]
check('4a sigma=(r13 r14)(r23 r24) permutes {+-g1..g5} => sigma(z) is a common zero of g1..g5 (NOT of g6)',
      all(any(h==g or h==-g for g in G[:5]) for h in SG[:5]))
CLOSED=[1,2,4,5,6,7,8,9,10,11,13,14]   # E_i built from g1..g5 only (E3, E12 use g6 via L6)
cs=flint.fmpz_mpoly_ctx.get(('a','s','p'),'lex'); A,S,P=cs.gens(); zs=0*A
from math import comb
def tosp(f):  # symmetric f in Z[a][r23,r24] -> f in Z[a][s,p], s=r23+r24, p=r23 r24
    f={(e[0],e[5],e[6]):c for e,c in f.to_dict().items()}; assert all(e[2]==e[3]==e[4]==0 for e in []) ; out={}
    while f:
        k=max(f,key=lambda k:(k[1],k[2],k[0])); ea,i,j=k; c=f[k]; assert i>=j
        m=i-j; out[(ea,m,j)]=out.get((ea,m,j),0)+c
        for tt in range(m+1):
            kk=(ea,j+m-tt,j+tt); v=f.get(kk,0)-c*comb(m,tt)
            if v: f[kk]=v
            else: f.pop(kk,None)
    return cs.from_dict(out)
def back(F):  # substitute s=r23+r24, p=r23 r24
    return sum((c*a**e[0]*(r23+r24)**e[1]*(r23*r24)**e[2] for e,c in F.to_dict().items()),zero)
dlt=r23-r24
SP={}
for i in CLOSED:
    e=E[i]; es=e+sw(e); ea_,rem=divmod(e-sw(e),dlt)
    assert rem.is_zero()
    SP['%ds'%i]=tosp(es); SP['%da'%i]=tosp(ea_)
    assert back(SP['%ds'%i])==es and back(SP['%da'%i])==ea_
for i in (3,12):   # g6-derived: only E_i(z)=0 is known, so use the symmetric norm E_i*sw(E_i)
    SP['%dn'%i]=tosp(E[i]*sw(E[i])); assert back(SP['%dn'%i])==E[i]*sw(E[i])
check('4b E_i+sw(E_i) and (E_i-sw(E_i))/(r23-r24) vanish at z; rewritten exactly in (s,p) (back-substitution checked); norms N3,N12 of the g6-derived E3,E12',True)
FF={'F1':A*P**2-S**2+2*P,'F2':A*P**2+S**2-2*P,'F3':A*P**3+S**3-3*S*P,'F4':A**2*P**4-S**4+4*S**2*P}
def divout(f,qs):
    for q in qs:
        while True:
            qq,r=divmod(f,q)
            if not r.is_zero() or qq.is_constant(): break
            f=qq
    return f
Hl={k:divout(v,[P,A]+list(FF.values())) for k,v in SP.items()}   # leaf: p,a,F1..F4 nonzero at z
Hb={k:divout(v,[P,A]) for k,v in SP.items()}                     # F-branches: only p,a divided
# T = Xn*Yn*Sn (= D^3 r13^2 r14^2 r34^2 at z, nonzero)
check('4c Xn,Yn,Sn in Z[a][r23,r24] and Xn*Yn*Sn symmetric',all(all(e[2]==e[3]==e[4]==0 for e in q.to_dict()) for q in (Xn,Yn,Sn))
      and sw(Xn*Yn*Sn)==Xn*Yn*Sn)
T=tosp(Xn*Yn*Sn); check('4c2 T(s,p) back-substitutes to Xn*Yn*Sn',back(T)==Xn*Yn*Sn)

Pst=(A**2-1)*P-3      # spurious value p0=3/(a^2-1)
Lsp=(A**2-1)*S+3*A    # spurious s0=-3a/(a^2-1)
def nonconst_factors(g,var):   # factors of g involving var (a-content dropped)
    return [(f,m) for f,m in g.factor()[1] if f.degrees()[var]>0]
def spurious(polys,name):
    """common zeros with P*(p0)=0: s0 is a root of gcd_s Res_p(f,P*) = L^k; and L | Res_p(T,P*) => T(s0,p0)=0, contradiction"""
    g=None
    for f in polys:
        R=f.resultant(Pst,'p'); g=R if g is None else g.gcd(R)
    fs=nonconst_factors(g,1)
    q,r=divmod(T.resultant(Pst,'p'),Lsp)
    check('%s: on P*=0 the common s-roots satisfy L=(a^2-1)s+3a=0, and L | Res_p(T,P*) (T=0: r13 r14 r34=0)'%name,
          len(fs)>0 and all(f==Lsp or f==-Lsp for f,m in fs) and r.is_zero())
spurious([Hl['10s'],Hl['10a']],'4d leaf spurious point')
# leaf: R1=Res_s(H10s,H10a) exact; R1' = R1 without a-content, p, P*; Res_p(R1',Res_s(H10s,H12n)) != 0 by specialisation
R1=Hl['10s'].resultant(Hl['10a'],'s')
R1p=cs.from_dict({(0,0,0):1})
for f,m in nonconst_factors(R1,2):
    if f!=P and f!=Pst and f!=-Pst: R1p*=f**m
def spec_coprime(f,g1,g2,a0,q):
    """True => Res_p(f, Res_s(g1,g2)) is a nonzero element of Z[a] (f in Z[a][p]; s-degrees and deg_p f preserved at a0 mod q)"""
    cq=flint.nmod_mpoly_ctx.get(('s','p'),ordering='lex',modulus=q)
    def sp(F):
        d={}
        for e,c in F.to_dict().items(): k=(e[1],e[2]); d[k]=(d.get(k,0)+int(c)*pow(a0,e[0],q))%q
        return cq.from_dict({k:v for k,v in d.items() if v})
    G1,G2,FF_=sp(g1),sp(g2),sp(f)
    if G1.degrees()[0]!=g1.degrees()[1] or G2.degrees()[0]!=g2.degrees()[1]: return False
    if FF_.degrees()[1]!=f.degrees()[2]: return False
    def up(F):
        d=F.to_dict(); cl=[0]*(max(k[1] for k in d)+1)
        for k,c in d.items(): cl[k[1]]=int(c)
        return flint.nmod_poly(cl,q)
    R2=G1.resultant(G2,'s')
    if R2.is_zero(): return False
    return up(FF_).gcd(up(R2)).degree()==0
check('4e leaf: Res_p(R1\', Res_s(H10s,H12n)) != 0 (specialisation a0=7 mod 1000003) => p0 in {0, 3/(a^2-1)}',
      spec_coprime(R1p,Hl['10s'],Hl['12n'],7,1000003))
# F-branches
for k in ('F1','F2','F3','F4'):
    g=None
    for j in ('1s','1a','3n','12n'):
        R=FF[k].resultant(Hb[j],'s'); g=R if g is None else g.gcd(R)
    fs=nonconst_factors(g,2)
    if k!='F3': check('4f branch %s=0: gcd_j Res_s(%s,H_j) is a monomial a^i p^j => p=0, impossible'%(k,k),all(f==P for f,m in fs))
    else:
        O=A**2*P**3-4
        check('4f branch F3=0: gcd factors only p, P*, a^2p^3-4',all(f in (P,Pst,-Pst,O,-O) for f,m in fs))
spurious([FF['F3'],Hb['1s'],Hb['1a']],'4g F3 spurious point')
# F3 & a^2p^3=4 => r23^3=r24^3=-2/a, r23!=r24 => r24^2+r23 r24+r23^2=0; E1 is nonzero at all such points:
c2=flint.fmpz_mpoly_ctx.get(('a','x','y'),'lex'); aa,X_,Y_=c2.gens()
f1=sum((c*aa**e[0]*X_**e[5]*Y_**e[6] for e,c in E[1].to_dict().items()),0*aa)
Nrm=f1.resultant(Y_**2+X_*Y_+X_**2,'y').resultant(aa*X_**3+2,'x')
check('4h F3 & a^2p^3=4: norm Res_x(a x^3+2, Res_y(y^2+xy+x^2, E1)) is a nonzero element of Z[a]',
      (not Nrm.is_zero()) and all(e[1]==e[2]==0 for e in Nrm.monoms()))
print('ALL PASS: no torus zero for transcendental a (b=1); cells 334 and 1767 excluded for generic masses' if all(OK) else 'FAIL')
sys.exit(0 if all(OK) else 1)
