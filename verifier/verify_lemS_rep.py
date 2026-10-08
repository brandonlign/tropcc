#!/usr/bin/env python3
"""Complete condition (ii) of Lemma S for target cells 81, 121 and 206: the 37 forms at (a,b)=(3,7) mod p=32003 have no common
zero in the torus over F_p-bar. Exact replay (python-flint, no Singular), plus nine independently checked B0 cofactor identities per target in data/lemS/b0/.
Steps: (0) forms from sys.json+fan.out (verify_lemS.forms), torus slice r13=1 (homogeneity), even vars r23,r24,r25 -> s=r^2;
 (1) g4: s25=-b(s23+s24); (2) g0: s23*c23+s24*c24=0 -> branch B0 (c23=c24=0 <=> x^3=y^3=1) or s=(lam P, lam Q), lam!=0, Z!=0;
 (3) lambda-free forms LF; branches phi1,phi2,phi3; (4) main: S^2,T^2,W^2 from q's, 1/T,1/W affine in 1/S (g9,g13), 1/S rational;
 (5) leaf: ES,EW,E27,E28 have no common zero outside the excluded set (resultant + extension-field gcds)."""
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
def mnorm(f):  # divide by gcd monomial and leading coefficient (monomials and constants are units on the torus)
    d=f.to_dict(); m=[min(e[i] for e in d) for i in range(len(VV))]
    g=C.from_dict({tuple(z-m[i] for i,z in enumerate(e)):c for e,c in d.items()})
    return g*pow(int(g.leading_coefficient()),-1,P0)

# The exceptional branch is mandatory: all nine finite-field identities are replayed.
assert '--noBR' not in sys.argv, 'partial checks cannot certify a cell'
import b0_verify,types
check('2c all nine B0 branches',b0_verify.verify(types.SimpleNamespace(**globals())))

# main: s23=lam P, s24=lam Q, lam != 0 (s23 != 0), Z != 0 by 2b
LF={}
for i,g in sorted(G1.items()):
    h=sub(g,s23=lam*P,s24=lam*Q)
    if h.is_zero() or deg(h,'r12')>0: continue
    ks=[k for k in range(deg(h,'lam')+1) if not coef(h,'lam',k).is_zero()]
    if len(ks)==1: LF[i]=mnorm(coef(h,'lam',ks[0]))   # lam^k and monomials are units
check('3a lambda-homogeneous forms -> lambda-free LF: %s'%sorted(LF),all(k in LF for k in list(K['sq'].values())+[K['DST'],K['DSW']]+list(K['leaf'])))
# ---------- 4. main branch: phi1 phi2 phi3 != 0 ----------
class Fq:   # fraction with gcd reduction (reduced denominators divide the unreduced ones, which are nonzero on the branch)
    def __init__(s,n,d=None):
        d=C.constant(1) if d is None else d
        assert not d.is_zero(), 'zero polynomial denominator'
        g=n.gcd(d) if not n.is_zero() else d
        s.n,s.d=(n//g,d//g) if not g.is_zero() else (n,d)
    def __add__(s,o): return Fq(s.n*o.d+o.n*s.d,s.d*o.d)
    def __sub__(s,o): return Fq(s.n*o.d-o.n*s.d,s.d*o.d)
    def __mul__(s,o): return Fq(s.n*o.n,s.d*o.d)
    def __truediv__(s,o): return Fq(s.n*o.d,s.d*o.n)
one=Fq(C.constant(1))
import leafchk
def xy(f):
    assert all(all(not e[j] for j in range(2,len(VV))) for e in f.to_dict()), 'non-bivariate leaf input'
    return leafchk.to_xy(f,P0,0,1)[0]
UNIT=x*y*Zx
def zfac(f):  # all irreducible factors of f (over F_p) divide x*y*Zx
    return all(UNIT%g==0 for g,e in f.factor()[1])
sig={}; c0s=[]
for v,k in K['sq'].items():
    f=LF[k]; c2,c0=coef(f,v,2),coef(f,v,0)
    others=[w for w in 'STW' if w!=v]
    good=f==c2*C.gens()[VV.index(v)]**2+c0 and all(deg(c2,w)<=0 and deg(c0,w)<=0 for w in others)
    if not zfac(c2):
        check('4a denominator of '+v+' cannot vanish on the torus branch',leafchk.leaf([xy(mnorm(c2)),xy(mnorm(c0))],[xy(f) for f in [x,y,P,Q,S25p]],P0))
        UNIT*=c2
    check('4a LF%d = c2 %s^2 + c0, c2 has verified nonzero factors, so %s^2 = -c0/c2'%(k,v,v),good and zfac(c2))
    sig[v]=Fq(-c0,c2); c0s.append(c0)
def tri(f,A,B):   # f = cAB*A^3B^3 + cA*A^3 + cB*B^3 ?
    VA,VB=C.gens()[VV.index(A)],C.gens()[VV.index(B)]
    cAB=coef(coef(f,A,3),B,3); cA=coef(coef(f,A,3),B,0); cB=coef(coef(f,A,0),B,3)
    return f==cAB*VA**3*VB**3+cA*VA**3+cB*VB**3,cAB,cA,cB
ok9,A9,k1,k2=tri(LF[K['DST']],'S','T'); ok13,A13,k3,k4=tri(LF[K['DSW']],'S','W')
print('  k1',k1,' k2',k2,' k3',k3,' k4',k4,' A9',A9,' A13',A13)
check('4b D_ST = A9 S^3T^3 + k1 S^3 + k2 T^3, A9 with only Z-factors',ok9 and zfac(A9))
check('4c D_SW = A13 S^3W^3 + k3 S^3 + k4 W^3, A13 with only Z-factors, k4 ~ k2',ok13 and zfac(A13) and prop(k4,k2))
phi1,phi2,phi3=k1,k2,k3
# divide by S^3T^3 (resp S^3W^3); 1/S^3=p/sigS etc. with p=1/S:  A9 + k1 q/sigT + k2 p/sigS = 0
F=lambda f: Fq(f)
al=Fq(C.constant(0))-F(A9)*sig['T']/F(k1); be=Fq(C.constant(0))-F(k2)*sig['T']/(sig['S']*F(k1))
ga=Fq(C.constant(0))-F(A13)*sig['W']/F(k3); de=Fq(C.constant(0))-F(k4)*sig['W']/(sig['S']*F(k3))
# q^2 sigT = 1 and p^2 = 1/sigS  =>  p = (1-(al^2+be^2/sigS) sigT)/(2 al be sigT)
p=(one-(al*al+be*be/sig['S'])*sig['T'])/(F(C.constant(2))*al*be*sig['T'])
q=al+be*p; r=ga+de*p
def num(fr): return fr.n
ES=num(p*p*sig['S']-one); EW=num(r*r*sig['W']-one); ET=num(q*q*sig['T']-one)
def ev(f):  # f(S=1/p,T=1/q,W=1/r) numerator
    dS,dT,dW=deg(f,'S'),deg(f,'T'),deg(f,'W'); out=C.from_dict({}); by={}
    for e,c in f.to_dict().items():
        k=(e[2],e[3],e[4]); m=list(e); m[2]=m[3]=m[4]=0; by.setdefault(k,{})[tuple(m)]=c
    for (i,j,k),t in by.items(): out+=C.from_dict(t)*p.d**i*p.n**(dS-i)*q.d**j*q.n**(dT-j)*r.d**k*r.n**(dW-k)
    return out
E27,E28=ev(LF[K['leaf'][0]]),ev(LF[K['leaf'][1]])
print('  leaf sizes',[ (len(e),deg(e,'x'),deg(e,'y')) for e in (ES,EW,E27,E28)],flush=True)
EXCL=[x,y,P,Q,S25p,phi1,phi2,phi3]+c0s+[p.n,p.d,q.n,q.d,r.n,r.d]
import leafchk
check('4d leaf: ES,EW,E27,E28 have no common zero in F_p-bar^2 off the excluded set',leafchk.leaf([xy(f) for f in (ES,EW,E27,E28)],[xy(f) for f in EXCL],P0))


# ---------- 5. branches phi_i = 0 by replay (alternative to lift certificates) ----------
cube=lambda fr: fr*fr*fr
def br_leaf(name,ph,Fs_fr,extra_ex):
    Fs=[ph]+[f.n for f in Fs_fr]
    ex=[x,y,P,Q,S25p]+c0s+extra_ex
    check('5 branch %s=0 (replay): no common zero off the excluded set'%name,leafchk.leaf([xy(f) for f in Fs],[xy(f) for f in ex],P0))
zero=Fq(C.constant(0))
# case logic: phi1=0: S^3=-k2/A9 (T!=0); if A13 S^3+k4=0 then k3 S^3=0 -> phi3=0: case phi1&phi3 below; else W^3 from LF13.
# phi3=0: S^3=-k4/A13 (W!=0); if A9 S^3+k2=0 then k1=0 -> case phi1&phi3; else T^3 from LF9.  phi2=0: k2=k4=0, T^3=-k1/A9, W^3=-k3/A13.
S3=zero-F(k2)/F(A9); W3=zero-F(k3)*S3/(F(A13)*S3+F(k4))
br_leaf('phi1',phi1,[S3*S3-cube(sig['S']),W3*W3-cube(sig['W'])],[A9,(F(A13)*S3+F(k4)).n,S3.n,W3.n])
# phi2=0 (k2=k4=0): T^3=-k1/A9, W^3=-k3/A13
T3=zero-F(k1)/F(A9); W3b=zero-F(k3)/F(A13)
br_leaf('phi2',phi2,[T3*T3-cube(sig['T']),W3b*W3b-cube(sig['W'])],[A9,A13,T3.n,W3b.n])
# phi3=0 (k3=0): S^3=-k4/A13 ; T^3 from LF9
S3c=zero-F(k4)/F(A13); T3c=zero-F(k1)*S3c/(F(A9)*S3c+F(k2))
br_leaf('phi3',phi3,[S3c*S3c-cube(sig['S']),T3c*T3c-cube(sig['T'])],[A13,(F(A9)*S3c+F(k2)).n,S3c.n,T3c.n])

# phi1=phi3=0: S^3=-k2/A9 and S^6=sigS^3
br_leaf('phi1&phi3',phi1,[Fq(phi3),S3*S3-cube(sig['S'])],[A9,S3.n])
print('verify_lemS_rep cell %d: ALL PASS'%CELL if all(OK) else 'FAIL')
