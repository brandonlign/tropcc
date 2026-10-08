# exact mod-p leaf check: no common zero of polys Fs in F_p-bar^2 outside V(prod EX); Fs, EX: nmod_mpoly in ctx (x,y)
import flint,time
def to_xy(f,p,ix,iy):
    assert all(all(not e[j] for j in range(len(e)) if j not in (ix,iy)) for e in f.to_dict()), 'unexpected extra variables'
    C2=flint.nmod_mpoly_ctx.get(('x','y'),p,'lex')
    return C2.from_dict({(e[ix],e[iy]):int(c) for e,c in f.to_dict().items()}),C2
def ycoeffs(f):  # dict xpower -> list of y-coeffs
    d={}
    for (i,j),c in f.to_dict().items():
        d.setdefault(i,{})[j]=int(c)
    return d
def leaf(Fs,EX,p,verbose=True):
    t0=time.time()
    r=Fs[0].resultant(Fs[1],'x')
    ry={e[1]:int(c) for e,c in r.to_dict().items()}; n=max(ry) if ry else -1
    R=flint.nmod_poly([ry.get(k,0) for k in range(n+1)],p)
    assert not R.is_zero(), 'resultant identically zero'
    M=flint.fmpz_mod_poly_ctx(p)
    fac=R.factor()[1]
    if verbose: print('Res deg',R.degree(),'distinct factors',len(fac),'max deg',max(f.degree() for f,e in fac),round(time.time()-t0,1),flush=True)
    YC=[ycoeffs(F) for F in Fs]; XC=[ycoeffs(E) for E in EX]
    bad=[]
    for f,e in sorted(fac,key=lambda t:t[0].degree()):
        K=flint.fq_default_ctx(modulus=M([int(c) for c in f.coeffs()]))
        PR=flint.fq_default_poly_ctx(K)
        def sp(yc):
            n=max(yc); return PR([K([yc[i].get(j,0) for j in range(max(yc[i])+1)]) if i in yc else K(0) for i in range(n+1)])
        g=None
        for yc in YC:
            q=sp(yc)
            if q.is_zero(): continue
            g=q if g is None else g.gcd(q)
            if g.degree()==0: break
        if g is None: bad.append((f.degree(),'allzero')); continue
        if g.degree()>0:
            for xc in XC:
                q=sp(xc)
                if q.is_zero(): g=PR([K(1)]); break      # excluded polynomial vanishes identically on this fibre
                while g.degree()>0:
                    h=g.gcd(q)
                    if h.degree()==0: break
                    g=g//h
        if g.degree()>0: bad.append((f.degree(),g.degree()))
    if verbose: print('leaf done',round(time.time()-t0,1),'bad',bad,flush=True)
    return not bad
if __name__=='__main__':
    import pickle
    from el2 import C as C0
    import der2
    P0=32003
    L={k:C0.from_dict(d) for k,d in pickle.load(open('leaf_spec.pkl','rb')).items()}
    sp=lambda f: f.compose(*([C0.constant(3),C0.constant(7)]+list(C0.gens())[2:]))
    dd=pickle.load(open('cub.pkl','rb'))
    a,b,lam,x,y=C0.gens()[:5]; X,Y=x**3,y**3
    ex=[sp(C0.from_dict(c0)) for v,(c2,c0) in dd['q'].items()]+[sp(getattr(der2,k)) for k in ('pn','pd','qn','qd','rn','rd')]
    ex+=[x,y,x**3-1,y**3-1,x**3-y**3,sp(b*X*Y-b*Y+X-Y),X*Y-2*X+Y,sp(b*X*Y-b*Y+X*Y-X)]
    Fs=[to_xy(L[k],P0,3,4)[0] for k in ('ES','EW','E27','E28')]
    EX=[to_xy(e,P0,3,4)[0] for e in ex]
    print(leaf(Fs,EX,P0))
