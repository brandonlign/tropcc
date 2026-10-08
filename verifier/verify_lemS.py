#!/usr/bin/env python3
"""Complete Lemma S hypothesis check for cells 81, 121 and 206.

Reconstruct the target's initial forms and the entire star from the fan;
require all original/augmented forms to be constant and U-free, and every
coefficient to survive at (3,7) modulo 32003. Independently replay the target
cell (including B0), and check every proper star cone by a modular cofactor
identity or the explicit cluster replay. Only a complete pass certifies a
cell's generic exclusion. No Singular computation is used here.
usage: verify_lemS.py 81 121 206
"""
import sys,os,json
import flint
import verify as V
HERE=os.path.dirname(os.path.abspath(__file__)); CD=os.path.join(V.N5,'lemS','cert')
rv,F=V.load_sys(); rv,F=V.augmentH(rv,F); rays,cells=V.read_fan(os.path.join(V.N5,'fan.out')); iU=rv.index('U')
def forms(c):
    Vs=[rays[i]+[0]*(len(rv)-len(rays[i])) for i in cells[c]]; w=[sum(x) for x in zip(*Vs)]; out={}
    for i,f in enumerate(F):
        sup=list(f); Fi=set(sup)
        for v in Vs:
            M=max(V.dot(v,m) for m in sup); Fi&={m for m in sup if V.dot(v,m)==M}
        Mw=max(V.dot(w,m) for m in sup); Fw={m for m in sup if V.dot(w,m)==Mw}
        if Fi and Fi==Fw and all(m[iU]==0 for m in Fw): out[i]={m:f[m] for m in Fw}
    return out
def spec(cf,p,a0,b0): return sum(v*pow(a0,ea,p)*pow(b0,eb,p) for (ea,eb),v in cf.items())%p
def check_cone(c,p0=None):
    fn=os.path.join(CD,'%d.cof'%c); L=open(fn).read().split('\n')
    assert L[0]=='cone %d'%c
    t=L[1].split(); p,a0,b0=int(t[1]),int(t[3]),int(t[5])
    if p0: assert (p,a0,b0)==p0, 'spec mismatch'
    vs=L[2][5:].split(','); assert vs[-1]=='u'; fl=[int(z) for z in L[3][6:].split(',')]
    fs=forms(c); assert set(fl)<=set(fs), 'cert uses a form that is not constant/U-free on the cone'
    ctx=flint.nmod_mpoly_ctx.get(tuple(vs),p,'lex'); n=len(vs); pos=[rv.index(v) for v in vs[:-1]]
    assert all(all(m[j]==0 for m in g for j in range(len(rv)) if j not in pos) for i,g in fs.items() if i in fl), 'form uses a variable outside vars'
    G=[]
    for i in fl:
        g=fs[i]; mins=[min(m[j] for m in g) for j in range(len(rv))]; d={}
        for m,cf in g.items():
            e=tuple(m[j]-mins[j] for j in pos)+(0,); d[e]=(d.get(e,0)+spec(cf,p,a0,b0))%p
        G.append(ctx.from_dict({e:c_ for e,c_ in d.items() if c_}))
    u=ctx.gens()[-1]; pr=ctx.constant(1)
    for g_ in ctx.gens()[:-1]: pr*=g_
    G.append(1-u*pr)
    T=[]; cur=None
    for line in L[4:]:
        if not line: continue
        if line[0]=='#': cur={}; T.append(cur); continue
        cs,es=line.split('|'); e=tuple(int(z) for z in es.split(',')); cur[e]=int(cs)%p
    assert len(T)==len(G)
    S=ctx.from_dict({})
    for t_,g_ in zip(T,G): S+=ctx.from_dict(t_)*g_
    return S==ctx.constant(1),(p,a0,b0)
def main(k):
    import subprocess,hashlib
    assert k in (81,121,206), 'structural replay supports only these three target cells'
    p0=(32003,3,7)
    fs=forms(k)
    # This strong check also guarantees every neighbouring certificate uses
    # only generators in Phi_C, as required by the statement of Lemma S.
    assert set(fs)==set(range(len(F))), 'not all original and augmented forms are available'
    bad=[(i,m) for i,g in fs.items() for m,cf in g.items() if spec(cf,*p0)==0]
    if bad:
        print('FAIL vanished coefficients',bad[:3]); return False
    st=[j for j,c in enumerate(cells) if set(cells[k])<=set(c)]
    print('cell',k,'all',len(fs),'coefficient supports survive: PASS',flush=True)
    ld=os.path.join(V.N5,'lemS','verification'); os.makedirs(ld,exist_ok=True)
    def replay(script,c):
        fn=os.path.join(ld,'%s_%d.log'%(script.removesuffix('.py'),c))
        with open(fn,'w') as stream:
            result=subprocess.run([sys.executable,os.path.join(HERE,script),str(c)],stdout=stream,stderr=subprocess.STDOUT)
        print('replay',script,c,'PASS' if result.returncode==0 else 'FAIL',flush=True)
        return result.returncode==0
    allok=replay('verify_lemS_rep.py',k)
    for c in st:
        if c==k: continue
        if c in (2174,2885,3494):
            ok=replay('cluster_mod.py',c)
        else:
            try:
                fn=os.path.join(CD,'%d.cof'%c)
                fl=[int(z) for z in open(fn).read().splitlines()[3][6:].split(',')]
                assert set(fl)<=set(fs), 'certificate uses a generator outside Phi_C'
                ok,_=check_cone(c,p0)
            except (AssertionError,FileNotFoundError,ValueError) as e:
                ok=False; print('cone',c,'ERROR',e,flush=True)
        allok=bool(allok and ok)
        if not ok: print('cone',c,'FAIL',flush=True)
    print('star(%d): %d cones; Lemma S hypotheses:'%(k,len(st)),'ALL PASS' if allok else 'FAIL',flush=True)
    # The marker is only a convenience for exploratory --fast coverage.
    marker=os.path.join(V.N5,'lemS','%d.ok'%k)
    if allok: open(marker,'w').write('all coefficient, cell, and star checks passed\n')
    elif os.path.exists(marker): os.unlink(marker)
    return allok
if __name__=='__main__':
    results=[main(int(k)) for k in sys.argv[1:]]
    sys.exit(0 if results and all(results) else 1)
