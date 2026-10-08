#!/usr/bin/env python3
"""Exact finite-field cluster replay for star cones 2174, 2885 and 3494.

All inputs are reconstructed at (a,b)=(3,7), p=32003. The Cramer and
division-free elimination steps yield fourteen bivariate equations, whose
common roots are excluded by an exact resultant/fiber-gcd leaf check.
This is a modular unit proof; generic exclusion requires the complete
star argument in verify_lemS.py. It is not a generic certificate by itself.
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
V=('a','b','r13','r14','r34','r23','r24'); ctx=flint.nmod_mpoly_ctx.get(V,32003,'lex')
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
TARGET=int(sys.argv[1]) if len(sys.argv)>1 and sys.argv[1].isdigit() else 2174
assert TARGET in (2174,2885,3494)
if TARGET==2174:
    raw=getforms(TARGET)
else:
    import transport,tempfile
    base={2885:2171,3494:2175}[TARGET]
    with tempfile.TemporaryDirectory(prefix='tropcc_cluster_') as td:
        fin=os.path.join(td,'input.cof'); fout=os.path.join(td,'output.cof')
        hdr=['cell %d'%base,'system F+G+H','divmono','subset 2,3,8,10,11,18','vars r14,r15,r34,r35,r45,u']
        open(fin,'w').write('\n'.join(hdr+['rhs prod^1','#1'])+'\n')
        transport.run(fin,TARGET,fout)
        lines=[l for l in open(fout).read().splitlines() if l and not l.startswith(('rhs ','#'))]
        result=VF.verify(None,lines=lines+[''],ret_gens=True)
        assert len(result)==4 and len(result[3])==7
        raw=[]
        for g in result[3][:-1]:
            d={}
            for e,v in g.to_dict().items():
                v=flint.fmpq(v); assert v.q==1 and e[7]==0
                d[(e[0],e[1],e[2],e[3],e[6],e[4],e[5])]=int(v.p)
            raw.append(d)
check('0 initial forms reconstructed for modular cluster cell %d'%TARGET,len(raw)==6)
G=[]
for d in raw:
    t={}
    for e,c in d.items():
        k=(0,0)+e[2:]; t[k]=(t.get(k,0)+c*pow(3,e[0],32003)*pow(7,e[1],32003))%32003
    G.append(ctx.from_dict({e:c for e,c in t.items() if c}))
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



Cxy=flint.nmod_mpoly_ctx.get(('x','y'),32003,'lex'); x,y=Cxy.gens()
def xy(f):
    assert all(all(e[j]==0 for j in range(5)) for e in f.to_dict())
    return Cxy.from_dict({(e[5],e[6]):c for e,c in f.to_dict().items()})
EE=[xy(f) for f in E.values()]
EX=[x,y,xy(Xn),xy(Yn),xy(Sn)]
def verify():
    import leafchk
    check('3 modular elimination leaf: all roots force a zero torus coordinate',leafchk.leaf(EE,EX,32003))
    return all(OK)
if __name__=='__main__':
    verify()
