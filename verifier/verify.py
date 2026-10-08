#!/usr/bin/env python3
"""Independent verifier for tropcc cell-exclusion certificates (no Singular, no sympy).
For certificate file C (cert/<k>.cof, produced by data/exact/certgen.py):
  1. read fan.out itself; cell k = k-th nonempty cone; rays v_1..v_s; w = sum v_j.
  2. read sys.json; set c=1 (all AC polys are homogeneous of degree 1 in (a,b,c), CM of degree 0).
  3. for each generator i in the subset: F_i = intersection_j argmax_{m in supp f_i} <v_j,m>; require |F_i|>=1 and
     F_i == argmax <w,m>  (then in_{w'} f_i = sum_{m in F_i} for every w' in relint(cell): CELL_LEMMA).
  4. g_i = in_w f_i with r_j := 1 for r_j not in the certificate's variable list ("fixed" coords).
     Check: homogeneity space H of {in_w f_i : i in subset} projects ONTO the fixed coordinates (rank = #fixed);
     then a torus point of V(in_w f_i) can be moved by H to have fixed coords = 1, so V(g) empty => V(in_w) empty.
  4b. optional header 'even v1,..': listed keep-vars occur only with even exponents in g_i; exponents are halved (s=r^2).
  5. check the polynomial identity  sum_i T_i*g_i + T_last*(1-u*prod keep) == 1  exactly in Q(a,b)[keep,u]
     (after clearing the lcm L of denominators: sum (L T_i) g_i == L, with L != 0 in Q[a,b]).
  => <in_w f_i> is the unit ideal of K[r^{+-1}], K = alg. closure of Q(a,b): no relint point of the cell is in Trop(X).
     Moreover the identity specialises to every (a,b) with L(a,b)!=0 (and the used initial coefficients not all zero).
"""
import sys,json,os
from fractions import Fraction
from math import gcd
import flint
HERE=os.path.dirname(os.path.abspath(__file__))
N5=os.path.join(HERE,'..','data')
def read_fan(fn):
    L=open(fn).read().split('\n'); i=L.index('RAYS'); rays=[]
    j=i+1
    while L[j].strip():
        rays.append([int(x) for x in L[j].split('#')[0].split()]); j+=1
    i=L.index('CONES'); cones=[]; j=i+1
    while j<len(L) and L[j].strip():
        t=L[j].split('#')[0].strip(); assert t[0]=='{' and t[-1]=='}',t
        cones.append([int(x) for x in t[1:-1].split()]); j+=1
    return rays,[c for c in cones if c]

import re
def parse_poly(s, ctx, names):
    """exact parser: sums of [rational coeff *] var^e products, optional outer parentheses."""
    s=s.replace(' ','').replace('**','^')
    while s.startswith('(') and s.endswith(')') and s.count('(')==1: s=s[1:-1]
    assert '(' not in s, s
    idx={v:i for i,v in enumerate(names)}; n=len(names); d={}
    if s[0] not in '+-': s='+'+s
    for sign,body in re.findall(r'([+-])([^+-]+)',s):
        c=Fraction(1); e=[0]*n
        for f in body.split('*'):
            if f[0].isdigit(): c*=Fraction(f); continue
            v,p=(f.split('^')+['1'])[:2]; e[idx[v]]+=int(p)
        if sign=='-': c=-c
        t=tuple(e); d[t]=d.get(t,0)+c
    return ctx.from_dict({k:flint.fmpq(v.numerator,v.denominator) for k,v in d.items() if v})
_SYS=None
def load_sys():
    global _SYS
    if _SYS is None:
        D=json.load(open(os.path.join(N5,'sys.json'))); rv=D['rvars']
        names=['a','b','c']+rv; F=[]
        for f in D['polys']:
            s=f.replace(' ','').replace('**','^')
            if s[0] not in '+-': s='+'+s
            terms={}   # r-exponent -> {(ea,eb): int}   with c=1
            for sign,body in re.findall(r'([+-])([^+-]+)',s):
                c=1; e=[0]*len(names)
                for x in body.split('*'):
                    if x[0].isdigit(): c*=int(x); continue
                    v,p=(x.split('^')+['1'])[:2]; e[names.index(v)]+=int(p)
                if sign=='-': c=-c
                re_=tuple(e[3:]); ab=(e[0],e[1])
                dd=terms.setdefault(re_,{}); dd[ab]=dd.get(ab,0)+c
            terms={m:{k:v for k,v in cf.items() if v} for m,cf in terms.items()}
            assert all(terms.values()), 'coefficient vanishes at c=1'
            F.append(terms)
        _SYS=(rv,F)
    return _SYS
def augment(rv,F):
    """Lemma U augmentation: coordinates (r_ij, U); F_i unchanged (U-exponent 0); f_35 = G = sum_{i<j} m_i m_j r_ij^2 - U,
    masses (a,a,b,b,c=1), pairs in the order of rvars (r12,r13,...,r45). U has weight 0 on every cell."""
    assert rv==['r12','r13','r14','r15','r23','r24','r25','r34','r35','r45']
    mab=[(1,0),(1,0),(0,1),(0,1),(0,0)]; pr=[(i,j) for i in range(5) for j in range(i+1,5)]
    FA=[{m+(0,):cf for m,cf in f.items()} for f in F]; g={}
    for t,(i,j) in enumerate(pr):
        e=[0]*11; e[t]=2; g[tuple(e)]={(mab[i][0]+mab[j][0],mab[i][1]+mab[j][1]):1}
    g[tuple([0]*10+[1])]={(0,0):-1}; FA.append(g)
    return rv+['U'],FA
def augmentH(rv,F):
    """S14: Lemma U + definition of U.  After augment(): f_36 = H = sum_{i<j} m_i m_j prod_{pairs != ij} r - U*prod_{all pairs} r
    (= (sum m_i m_j / r_ij - U) * prod r_ij, which vanishes on every real CC since U is the potential)."""
    rv,F=augment(rv,F)
    mab=[(1,0),(1,0),(0,1),(0,1),(0,0)]; pr=[(i,j) for i in range(5) for j in range(i+1,5)]; h={}
    for t,(i,j) in enumerate(pr):
        e=[1]*10+[0]; e[t]=0; h[tuple(e)]={(mab[i][0]+mab[j][0],mab[i][1]+mab[j][1]):1}
    h[tuple([1]*11)]={(0,0):-1}; F=F+[h]
    return rv,F
def dot(u,v): return sum(x*y for x,y in zip(u,v))

def nullspace_Q(rows,n):
    """rational nullspace basis of the matrix with given integer rows (n columns)."""
    if not rows: return [[1 if i==j else 0 for i in range(n)] for j in range(n)]
    M=flint.fmpq_mat(len(rows),n,[x for r in rows for x in r])
    R,rank=M.rref(); piv=[]; r=0
    for c in range(n):
        if r<rank and R[r,c]!=0: piv.append(c); r+=1
    free=[c for c in range(n) if c not in piv]; B=[]
    for f in free:
        v=[flint.fmpq(0)]*n; v[f]=flint.fmpq(1)
        for i,p in enumerate(piv): v[p]=-R[i,f]
        B.append(v)
    return B
def verify(fn, verbose=False, lines=None, ret_gens=False):
    """ret_gens=True (verify_split.py): after the support/torus checks return (k,ctx,names,G) instead of reading cofactors"""
    rv,F=load_sys(); rays,cells=read_fan(os.path.join(N5,'fan.out'))
    L=lines if lines is not None else open(fn).read().split('\n')
    aug = L[1].strip() in ('system F+G','system F+G+H')
    if aug: rv,F=(augmentH if L[1].strip()=='system F+G+H' else augment)(rv,F); L=L[:1]+L[2:]
    divm = L[1].strip()=='divmono'     # each g_i divided by the gcd monomial of its terms (a unit of the Laurent ring)
    if divm: L=L[:1]+L[2:]
    k=int(L[0].split()[1]); S=[int(x) for x in L[1].split()[1].split(',')]
    vars_=L[2].split()[1].split(','); assert len(set(vars_))==len(vars_)
    even=[]
    if L[3].startswith('even '): even=L[3][5:].strip().split(','); L=L[:3]+L[4:]   # s_v=r_v^2 substitution
    rhs=None
    if L[3].startswith('rhs '): rhs=L[3][4:].strip(); L=L[:3]+L[4:]
    rab = vars_[-1]=='u'
    assert rab or (rhs is not None and rhs.startswith('prod^')), 'need Rabinowitsch u or rhs prod^N'
    keep=vars_[:-1] if rab else vars_; assert all(v in rv for v in keep)
    fixj=[j for j,v in enumerate(rv) if v not in keep]
    V=[rays[i]+[0]*(len(rv)-len(rays[i])) for i in cells[k]]; w=[sum(x) for x in zip(*V)]
    # 3. cell-lemma support check
    forms=[]
    for i in S:
        sup=list(F[i]); Fi=set(sup)
        for v in V:
            M=max(dot(v,m) for m in sup); Fi&={m for m in sup if dot(v,m)==M}
        Mw=max(dot(w,m) for m in sup); Fw={m for m in sup if dot(w,m)==Mw}
        if not Fi or Fi!=Fw: return k,False,'support not constant on cell for f%d'%i
        forms.append(sorted(Fw))
    # 4. torus slice: projection of homogeneity space onto fixed coords must be onto
    rows=[[x-y for x,y in zip(f[0],m)] for f in forms for m in f[1:]]
    H=nullspace_Q(rows,len(rv))
    if fixj:
        P=flint.fmpq_mat(len(H),len(fixj),[h[j] for h in H for j in fixj]) if H else None
        if P is None or P.rank()!=len(fixj): return k,False,'torus slice invalid'
    # 5. exact identity  sum_i (L*T_i)*g_i == L  in Q[a,b][keep,u]
    names=['a','b']+vars_; ctx=flint.fmpq_mpoly_ctx.get(tuple(names),'lex'); nv=len(names)
    kpos=[rv.index(v) for v in keep]
    G=[]
    for Fw,i in zip(forms,S):
        d={}
        for m in Fw:
            for (ea,eb),c in F[i][m].items():
                e=(ea,eb)+tuple(m[j] for j in kpos)+((0,) if rab else ()); d[e]=d.get(e,0)+c
        if divm:
            mn=[min(e[2+t] for e in d) for t in range(len(kpos))]
            d={e[:2]+tuple(e[2+t]-mn[t] for t in range(len(kpos)))+e[2+len(kpos):]:c for e,c in d.items()}
        if even:   # every exponent of an 'even' var must be even; halve it (g(r)=g'(r^2): a unit identity for g' in s pulls back to one for g)
            ep=[2+keep.index(v) for v in even]
            assert all(e[j]%2==0 for e in d for j in ep), 'even var with odd exponent'
            d={tuple(x//2 if j in ep else x for j,x in enumerate(e)):c for e,c in d.items()}
        G.append(ctx.from_dict({e:c for e,c in d.items() if c}))
    ones=[0]*nv
    one=ctx.from_dict({tuple(ones):1})
    if rab: G.append(one-ctx.from_dict({tuple([0,0]+[1]*(nv-2)):1}))
    if rhs is None: R0=one
    elif rhs.startswith('prod^'): R0=ctx.from_dict({tuple([0,0]+[int(rhs[5:])]*(nv-2)):1})
    else: R0=parse_poly(rhs,ctx,names); assert all(e[2:]==(0,)*(nv-2) for e in R0.to_dict()) and not R0.is_zero()
    if ret_gens: return k,ctx,names,G
    # read cofactors
    T=[[] for _ in G]; cur=None; cache={}
    def P(s):
        if s not in cache: cache[s]=parse_poly(s,ctx,names)
        return cache[s]
    for line in L[3:]:
        if not line.strip(): continue
        if line[0]=='#': cur=int(line[1:])-1; continue
        num,den,ex=line.split('|'); e=tuple([0,0]+[int(x) for x in ex.split(',')])
        assert len(e)==nv
        T[cur].append((num,den,e))
    if cur!=len(G)-1: return k,False,'cofactor count mismatch'
    if len(G)!=len(S)+rab: return k,False,'generator count mismatch'
    Lc=one
    for dens in {den for t in T for _,den,_ in t}:
        d=P(dens); assert not d.is_zero()
        assert all(e[2:]==(0,)*(nv-2) for e in d.to_dict()), 'denominator must lie in Q[a,b]'
        Lc=Lc*d/Lc.gcd(d) if not d.is_constant() else Lc
    tot=ctx.from_dict({})
    for Ti,g in zip(T,G):
        acc={}
        for num,den,e in Ti:
            q=P(num)*Lc/P(den)
            for ee,c in q.to_dict().items():
                t=tuple(x+y for x,y in zip(ee,e)); acc[t]=acc.get(t,0)+c
        tot+=ctx.from_dict({t:c for t,c in acc.items() if c})*g
    ok=(tot==Lc*R0)
    return k,ok,('identity OK%s; L=%s'%(' [F+G]' if aug else '',(Lc*R0).str()) if ok else 'IDENTITY FAILS')
if __name__=='__main__':
    bad=0
    for fn in sys.argv[1:]:
        try: k,ok,msg=verify(fn)
        except Exception as e: k,ok,msg=fn,False,'malformed certificate: %r'%e
        print(k,'PASS' if ok else 'FAIL',msg[:120],flush=True); bad+=not ok
    print('failures',bad); sys.exit(1 if bad else 0)
