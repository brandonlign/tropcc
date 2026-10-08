#!/usr/bin/env python3
"""S15: split-tree (case-split) certificates.  File <k>.split:
  5 header lines exactly as in a verify.py .cof with Rabinowitsch u (cell / system / divmono / subset / vars ...,u)
  then blocks
    node <path> split <p1> ; <p2>     (path '' written as 'R'; children <path>.1, <path>.2)
    node <path> leaf
  each node followed by 'rhs <L>' (L in Q[a,b]\\0) and cofactor sections '#1'..'#m' (lines num|den|exp, exp over vars incl u),
  m = #base generators (verify.py, incl. 1-u*prod) + #extra generators on the path (the p_j chosen at the ancestors).
  Identity checked exactly in Q[a,b][vars]:  sum_i T_i*g_i == L*p1*p2 (split) or == L (leaf).
Logic: K_path = <base gens, extras>.  leaf: L in K_path, L a nonzero constant over Q(a,b) => 1 in K_path.
  split: p1*p2 in K_path, 1 in K_path+(p1), 1 in K_path+(p2) => 1=(k1+c1 p1)(k2+c2 p2) in K_path.
Hence 1 in K_R = the ideal of verify.py, i.e. the same exclusion statement as a .cof certificate; masses must avoid zeros of all L.
The tree must be complete (every split has both children)."""
import sys,flint
if hasattr(sys,"set_int_max_str_digits"): sys.set_int_max_str_digits(0)
from verify import verify, parse_poly
def check(fn):
    L=open(fn).read().split('\n'); hdr=L[:5]
    R=verify(fn,lines=hdr+['#1'],ret_gens=True)
    if len(R)==3: return R          # support / torus-slice failure
    k,ctx,names,G=R
    nv=len(names); assert names[-1]=='u'
    nodes={}; cur=None; sec=None
    for line in L[5:]:
        line=line.strip()
        if not line: continue
        if line.startswith('node '):
            w=line.split(None,3); path=w[1]; kind=w[2]; assert path not in nodes
            ps=[parse_poly(x.strip(),ctx,names) for x in w[3].split(';')] if kind=='split' else []
            assert (kind=='split' and len(ps)==2) or (kind=='leaf' and not ps)
            cur=nodes[path]={'kind':kind,'ps':ps,'T':{},'rhs':None}; sec=None
        elif line.startswith('rhs '):
            r=parse_poly(line[4:].strip(),ctx,names); assert not r.is_zero() and all(e[2:]==(0,)*(nv-2) for e in r.to_dict()); cur['rhs']=r
        elif line[0]=='#': sec=int(line[1:])-1; cur['T'].setdefault(sec,ctx.from_dict({}))
        else:
            num,den,ex=line.split('|'); e=tuple([0,0]+[int(x) for x in ex.split(',')]); assert len(e)==nv
            d=parse_poly(den,ctx,names); assert all(t[2:]==(0,)*(nv-2) for t in d.to_dict())
            q=parse_poly(num,ctx,names)*ctx.from_dict({e:1})
            assert all(t[2:]==(0,)*(nv-2) for t in parse_poly(num,ctx,names).to_dict()), 'numerator must lie in Q[a,b]'
            if not d.is_constant(): return k,False,'non-constant denominator not supported (put it in rhs)'
            cur['T'][sec]+=q/d
    Ls=[]
    def walk(path,extras):
        if path not in nodes: return 'missing node '+path
        n=nodes[path]; gens=G+extras
        if n['rhs'] is None or any(i<0 or i>=len(gens) for i in n['T']): return 'bad node '+path
        tot=ctx.from_dict({})
        for i,t in n['T'].items(): tot+=t*gens[i]
        tgt=n['rhs']*(n['ps'][0]*n['ps'][1] if n['kind']=='split' else 1)
        if tot!=tgt: return 'identity fails at node '+path
        Ls.append(n['rhs'].str())
        if n['kind']=='split':
            for j in (0,1):
                r=walk(path+'.%d'%(j+1),extras+[n['ps'][j]])
                if r: return r
        return None
    r=walk('R',[])
    if r: return k,False,r
    if len(Ls)!=len(nodes): return k,False,'unreachable nodes'
    return k,True,'split tree OK (%d nodes); L=%s'%(len(nodes),' ; '.join(x[:40] for x in Ls))
if __name__=='__main__':
    bad=0
    for fn in sys.argv[1:]:
        try: k,ok,msg=check(fn)
        except Exception as e: k,ok,msg=fn,False,'malformed: %r'%e
        print(k,'PASS' if ok else 'FAIL',msg[:200],flush=True); bad+=not ok
    print('failures',bad); sys.exit(1 if bad else 0)
