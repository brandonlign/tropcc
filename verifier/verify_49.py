#!/usr/bin/env python3
"""Elimination certificate for the cluster 49-group (cells 49, 2171, 2174, 2175; S24). No Singular; exact (python-flint;
sympy only for one 4x4 subresultant). Claim: for generic (a,b) the 6 initial forms (subset 2,3,8,10,11,18, divmono, keep vars
r14,r15,r34,r35,r45) have no common zero on the complex torus. Pointwise proof over K = algebraic closure of Q(a,b)
(a,b algebraically independent), z a hypothetical common zero in (K^*)^5:
  0. forms via verify.py (support constancy + torus slice) for the four cells: identical.
  1-2. same division-free elimination as verify_334.py with roles (r13,r14,r34,r23,r24) <- (r14,r15,r45,r34,r35):
     E1..E14 in Z[a,b][x,y], x=r34, y=r35, vanishing at z.  Xn=-2x^3 Xq, Yn=-2b y^3 Yq (D r^2 = Xn etc) => Xq(z),Yq(z) != 0.
  3. leaf (all Cramer-Delta factors F0..F4 nonzero at z): H10,H13,H12 = cores of E10,E13,E12.  R1=Res_x(H10,H13) exact;
     Res_y(R1',Res_x(H10,H12)) != 0 by specialisation => phi(y)=0, phi=(ay+1)^3-b^2(ay^3+1) (irreducible).
  4. over L=Q(a,b)[y]/(phi): psi=-B/A from the first subresultant of Xq,Yq; Xq(psi)=0, (x-psi)^2 | H10, H13 and
     deg gcd_L(H10,H13) <= 2 (specialisation) => x(z)=psi => Xq(z)=0, contradiction.
  5. branches F_i(z)=0: gcds of Res_x(F_i, E_j); F0..F3 leave lc_x(F_i)=0 (then F_i = +-b y^2 != 0); F4 leaves phi, where
     A^3 F4 = (Ax+B) Q + r with r = 0 in L and Res_x(Q,E1) != 0 in L (specialisation) => x(z)=psi again.
"""
import sys,os,time,pickle
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
def getforms(k):
    L=['cell %d'%k,'system F+G+H','divmono','subset 2,3,8,10,11,18','vars r14,r15,r34,r35,r45,u','']
    r=VF.verify(None,lines=L,ret_gens=True)
    assert len(r)==4, r
    _,c0,names,G=r
    assert names==['a','b','r14','r15','r34','r35','r45','u']
    out=[]
    for g in G[:6]:
        d={}
        for e,v in g.to_dict().items():
            v=flint.fmpq(v); assert v.q==1 and e[7]==0
            d[(e[0],e[1],e[2],e[3],e[6],e[4],e[5])]=int(v.p)   # (r14,r15,r45,r34,r35) -> roles (r13,r14,r34,r23,r24)
        out.append(d)
    return out
FS={k:getforms(k) for k in (49,2171,2174,2175)}
check('0a forms of cells 49, 2171, 2174, 2175 (F+G+H, subset 2,3,8,10,11,18) coincide',all(FS[k]==FS[49] for k in FS))
G=[ctx.from_dict(d) for d in FS[49]]
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
check('2b L1=al r14+be, L2=ga r13+de, L6=al6 r14+be6 with coefficients in Z[a,b][r23,r24]',
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
check('2c M1=mu r14+nu, M2=m2 r13+n2 (coefficients in Z[a,b][r23,r24])',M1==mu*r14+nu and M2==m2*r13+n2
      and all(all(e[2]==e[3]==e[4]==0 for e in q.to_dict()) for q in (mu,nu,m2,n2)))
E[10]=strip(nu**2*D-mu**2*Yn); E[11]=strip(mu*be-al*nu); E[12]=strip(mu*be6-al6*nu); E[13]=strip(n2**2*D-m2**2*Xn); E[14]=strip(m2*de-ga*n2)
check('2d E1..E14 nonzero, in Z[a,b][r23,r24] (%.0fs)'%(time.time()-t0),
      all(not e.is_zero() and all(x[2]==x[3]==x[4]==0 for x in e.to_dict()) for e in E.values()))


# ---------- 3. leaf ----------
C4=flint.fmpz_mpoly_ctx.get(('a','b','x','y'),'lex'); A4,B4,x,y=C4.gens()
def to4(f): return sum((c*A4**e[0]*B4**e[1]*x**e[5]*y**e[6] for e,c in f.to_dict().items()),0*A4)
Xq=A4**2*x**4*y**4+A4*B4*x**3*y**4+A4*B4*x*y**4+B4**2*y**4-x**4
Yq=A4**2*x**4*y**4+A4*x**4*y**3+A4*x**4*y-B4**2*y**4+x**4
Sq=A4**2*x**4*y**4-B4**2*y**4-B4*x**3*y-B4*x*y**3-x**4
check('3a Xn=-2x^3 Xq, Yn=-2b y^3 Yq, Sn=2a x^3y^3 Sq (so Xq,Yq,Sq are nonzero at z: D r^2 = Xn etc, r != 0)',
      to4(Xn)==-2*x**3*Xq and to4(Yn)==-2*B4*y**3*Yq and to4(Sn)==2*A4*x**3*y**3*Sq)
FF=[A4*x**2*y**2-B4*y**2-x**2,A4*x**2*y**2-B4*y**2+x**2,A4*x**2*y**2+B4*y**2-x**2,A4*x**2*y**2+B4*y**2+x**2,A4*x**3*y**3+B4*y**3+x**3]
def divout(f,qs):
    for q in qs:
        while True:
            qq,r=divmod(f,q)
            if not r.is_zero() or qq.is_constant(): break
            f=qq
    return f
Hb={i:divout(to4(e)//to4(e).term_content(),[Xq,Yq,Sq]) for i,e in E.items()}   # Xq,Yq,Sq, monomials removed
H={i:divout(f,FF) for i,f in Hb.items()}                                       # leaf: F0..F4 removed too
phi=(A4*y+1)**3-B4**2*(A4*y**3+1)
c,fl=phi.factor()
check('3b phi=(ay+1)^3-b^2(ay^3+1) irreducible in Z[a,b,y] and primitive in y (=> irreducible over Q(a,b))',
      len(fl)==1 and fl[0][1]==1 and abs(int(c))==1)
t0=time.time()
R1=H[10].resultant(H[13],'x')
R1p=divout(R1,[y,phi])
check('3c R1=Res_x(H10,H13) = (content)*y^k*phi^m*R1\' computed exactly (%.0fs)'%(time.time()-t0),(not R1p.is_zero()))
def specialize(f,a0,b0,q,ctxq,idx):   # f in C4 -> nmod_mpoly in ctxq over the variables idx (subset of (x,y))
    d={}
    for e,c in f.to_dict().items():
        k=tuple(e[2+i] for i in idx); d[k]=(d.get(k,0)+int(c)*pow(a0,e[0],q)*pow(b0,e[1],q))%q
    return ctxq.from_dict({k:v for k,v in d.items() if v})
def upoly(F,i,q):   # nmod_mpoly -> nmod_poly in variable i
    d=F.to_dict(); cl=[0]*(max(k[i] for k in d)+1)
    for k,c in d.items(): cl[k[i]]=(cl[k[i]]+int(c))%q
    return flint.nmod_poly(cl,q)
A0,B0,QP=3,7,32003
cq=flint.nmod_mpoly_ctx.get(('x','y'),ordering='lex',modulus=QP)
h10,h12,r1p=[specialize(f,A0,B0,QP,cq,(0,1)) for f in (H[10],H[12],R1p)]
ok=h10.degrees()[0]==H[10].degrees()[2] and h12.degrees()[0]==H[12].degrees()[2] and r1p.degrees()[1]==R1p.degrees()[3]
R2s=h10.resultant(h12,'x')
check('3d Res_y(R1\', Res_x(H10,H12)) != 0: degrees preserved at (a,b)=(3,7) mod 32003 and gcd of specialisations = 1 => y(z) root of phi',
      ok and (not R2s.is_zero()) and upoly(r1p,1,QP).gcd(upoly(R2s,1,QP)).degree()==0)

# ---------- 4. the spurious points over L=Q(a,b)[y]/(phi) ----------
C3=flint.fmpz_mpoly_ctx.get(('a','b','y'),'lex'); a3,b3,y3=C3.gens(); Z3=0*a3
phi3=(a3*y3+1)**3-b3**2*(a3*y3**3+1); lcp=a3**3-a3*b3**2
class L:   # element num/lcp^k of L, num in Z[a,b][y] with deg_y<=2 (exact representation; zero test exact since lcp is a unit of L)
    def __init__(s,num,k=0):
        while not num.is_zero() and num.degrees()[2]>=3:
            d=num.degrees()[2]; top=C3.from_dict({(e[0],e[1],0):c for e,c in num.to_dict().items() if e[2]==d})
            num=num*lcp-top*y3**(d-3)*phi3; k+=1
        s.n=num; s.k=k
    def __add__(s,o): m=max(s.k,o.k); return L(s.n*lcp**(m-s.k)+o.n*lcp**(m-o.k),m)
    def __mul__(s,o): return L(s.n*o.n,s.k+o.k)
    def zero(s): return s.n.is_zero()
check('4a phi has leading y-coefficient lcp=a^3-ab^2',phi3.degrees()[2]==3 and C3.from_dict({(e[0],e[1],0):c for e,c in phi3.to_dict().items() if e[2]==3})==lcp)
def xco(f):   # C4 poly -> list of C3 coefficients of x^k
    co={}
    for e,c in f.to_dict().items(): co.setdefault(e[2],{}); kk=(e[0],e[1],e[3]); co[e[2]][kk]=co[e[2]].get(kk,0)+c
    return [C3.from_dict(co.get(k,{})) for k in range(max(co)+1)]
import sympy as sp
sa,sb,sx,sy=sp.symbols('a b x y')
def sym(f): return sp.sympify(str(f).replace('^','**'),locals=dict(a=sa,b=sb,x=sx,y=sy))
S1=[s for s in sp.subresultants(sp.Poly(sym(Xq),sx),sp.Poly(sym(Yq),sx)) if s.degree()==1][0]
def from_sym(e): P=sp.Poly(sp.expand(e),sa,sb,sy); return C3.from_dict({m:int(c) for m,c in P.terms()})
Acf,Bcf=[from_sym(c) for c in S1.all_coeffs()]     # S1 = Acf*x + Bcf, psi = -Bcf/Acf
check('4b A (lc of the first subresultant of Xq,Yq) is nonzero in L',not L(Acf).zero())
def hom(cs):   # A^n f(-B/A) in L, f = sum cs[k] x^k
    n=len(cs)-1; pa=[L(Z3+1)]; pb=[L(Z3+1)]
    for _ in range(n): pa.append(pa[-1]*L(Acf)); pb.append(pb[-1]*L(-Bcf))
    t=L(Z3)
    for k,c in enumerate(cs): t=t+L(c)*(pb[k]*pa[n-k])
    return t
def der(cs): return [k*c for k,c in enumerate(cs)][1:]
check('4c Xq(psi,y)=0 in L',hom(xco(Xq)).zero())
t0=time.time()
for j in (10,13):
    cs=xco(H[j]); cond=hom(cs).zero() and hom(der(cs)).zero()
    check('4d H%d(psi)=dH%d/dx(psi)=0 in L, i.e. (x-psi)^2 | H%d (%.0fs)'%(j,j,j,time.time()-t0),cond)
# upper bound deg gcd_L(H10,H13) <= 2: at a root y0 of phi mod q, lc_x preserved and gcd degree 2 (=> psc_2 not in (phi), Gauss)
rts=[int(-f.coeffs()[0])%QP for f,m in upoly(specialize(phi,A0,B0,QP,cq,(0,1)),1,QP).factor()[1] if f.degree()==1]
def at_y0(f,y0):
    cs=xco(f); return flint.nmod_poly([sum(int(c)*pow(A0,e[0],QP)*pow(B0,e[1],QP)*pow(y0,e[2],QP) for e,c in cc.to_dict().items())%QP for cc in cs],QP)
okb=False
for y0 in rts:
    u10,u13=at_y0(H[10],y0),at_y0(H[13],y0)
    if u10.degree()==H[10].degrees()[2] and u13.degree()==H[13].degrees()[2] and u10.gcd(u13).degree()==2: okb=True; break
check('4e deg gcd_L(H10,H13) <= 2 (specialisation at a root of phi mod 32003) => gcd_L = (x-psi)^2 => x(z)=psi, Xq(z)=0: contradiction',okb)

# ---------- 5. branches F_i(z)=0 (Cramer-Delta factors) ----------
BR=[1,2,3]   # E1,E2,E3 (Xq,Yq,Sq and monomials removed)
for i,F in enumerate(FF):
    g=None
    for j in BR:
        R=F.resultant(Hb[j],'x'); g=R if g is None else g.gcd(R)
    fs=[f for f,m in g.factor()[1] if f.degrees()[3]>0]
    if i<4:
        lc=sum((c*A4**e[0]*B4**e[1]*y**e[3] for e,c in F.to_dict().items() if e[2]==2),0*A4)
        c0=sum((c*A4**e[0]*B4**e[1]*y**e[3] for e,c in F.to_dict().items() if e[2]==0),0*A4)
        check('5a branch F%d=0: gcd_j Res_x(F%d,E_j) has only the factors y and lc_x(F%d); lc_x=0 forces F%d = %s != 0'%(i,i,i,i,c0),
              all(f==y or f==lc or f==-lc for f in fs) and F==lc*x**2+c0 and len(c0)==1)
    else:
        check('5b branch F4=0: gcd_j Res_x(F4,E_j) has only the factors y and phi',all(f==y or f==phi or f==-phi for f in fs))
        # pseudo-division A^3 F4 = (A x + B) Q + r  (A,B from step 4), r = hom(F4) must vanish in L
        cs=xco(F); n=len(cs)-1
        Ac4=sum((c*A4**e[0]*B4**e[1]*y**e[2] for e,c in Acf.to_dict().items()),0*A4); Bc4=sum((c*A4**e[0]*B4**e[1]*y**e[2] for e,c in Bcf.to_dict().items()),0*A4)
        # synthetic pseudo-division by (A x + B): Q = sum q_k x^k
        num=Ac4**n*F; Q=0*A4; cur=num
        for k in range(n-1,-1,-1):
            ck=sum((c*A4**e[0]*B4**e[1]*y**e[3] for e,c in cur.to_dict().items() if e[2]==k+1),0*A4)
            qk,rr=divmod(ck,Ac4); assert rr.is_zero()
            Q+=qk*x**k; cur=cur-qk*x**k*(Ac4*x+Bc4)
        r=cur; assert r.degrees()[2]==0
        rL=L(sum((c*a3**e[0]*b3**e[1]*y3**e[3] for e,c in r.to_dict().items()),Z3))
        check('5c A^3 F4 = (Ax+B) Q + r exactly, r = 0 in L',num==(Ac4*x+Bc4)*Q+r and rL.zero())
        okq=False
        for y0 in rts:
            uq,u1=at_y0(Q,y0),at_y0(Hb[1],y0)
            if uq.degree()==Q.degrees()[2] and u1.degree()==Hb[1].degrees()[2] and uq.gcd(u1).degree()==0: okq=True; break
        check('5d Res_x(Q,E1) nonzero in L (specialisation at a root of phi mod 32003) => x(z)=psi => Xq(z)=0: contradiction',okq)
print('ALL PASS: no torus zero for generic (a,b); cells 49, 2171, 2174, 2175 excluded' if all(OK) else 'FAIL')
sys.exit(0 if all(OK) else 1)
