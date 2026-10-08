"""Replay nine small cofactor identities over F_32003[t]/(t^2+t+1).

Inputs are regenerated from the fan/system and the checked torus reduction,
not trusted from Singular output. No Groebner computation is used here.
"""
import os
import flint

def verify(M, directory=None):
    p=M.P0
    from math import isqrt
    assert (M.P0,M.A0,M.B0)==(32003,3,7)
    assert p%3==2 and all(p%d for d in range(2,isqrt(p)+1))
    if M.CELL in (81,121):
        assert M.prop(M.Q,M.x**3*(M.y**3-1))
        assert M.prop(divmod(M.P,M.y**3-1)[1],M.x**3-1)
        rootx,rooty=1,1
    elif M.CELL==206:
        assert M.prop(M.Q,M.y**3*(M.x**3-1))
        assert M.prop(divmod(M.P,M.x**3-1)[1],13*M.y**3+1)
        rootx,rooty=1,pow(-pow(13,-1,p)%p,pow(3,-1,p-1),p)
    else:
        raise ValueError('unsupported cell')
    directory=directory or os.path.join(M.V.N5,'lemS','b0',str(M.CELL))
    names=('t','S','T','W','s23','s24','u')
    ctx=flint.nmod_mpoly_ctx.get(names,p,'lex')
    t,S,T,W,s23,s24,u=ctx.gens()
    modulus=t*t+t+1
    def reduce(f): return divmod(f,modulus)[1]
    def evaluate(f,ix,iy):
        d={}
        for e,c in f.to_dict().items():
            assert all(e[k]==0 for k in (7,8,9,10))
            v=int(c)*pow(rootx,e[0],p)*pow(rooty,e[1],p)%p
            key=((ix*e[0]+iy*e[1])%3,)+e[2:7]+(0,)
            d[key]=(d.get(key,0)+v)%p
        return reduce(ctx.from_dict(d))
    forms=[(i,g) for i,g in sorted(M.G1.items()) if M.deg(g,'r12')<=0 and not g.is_zero()]
    coefficient_ctx=flint.fmpq_mpoly_ctx.get(('t',),'lex')
    s25=-(M.l23*M.s23+M.l24*M.s24)*M.il25
    for ix in range(3):
        for iy in range(3):
            # Every root of the two cube equations is included in these nine cases.
            assert evaluate(M.P,ix,iy).is_zero() and evaluate(M.Q,ix,iy).is_zero()
            fn=os.path.join(directory,'%d_%d.cof'%(ix,iy))
            lines=open(fn).read().splitlines()
            assert lines[0]=='cell %d roots %d %d p %d'%(M.CELL,ix,iy,p)
            assert lines[1]=='vars '+','.join(names[1:])
            assert lines[2]=='forms '+','.join(str(i) for i,g in forms)
            polys=[evaluate(g,ix,iy) for i,g in forms]
            polys.append(1-u*S*T*W*s23*s24*evaluate(s25,ix,iy))
            cof=[]
            current=None
            for line in lines[3:]:
                if not line: continue
                if line.startswith('#'):
                    assert int(line[1:])==len(cof)+1
                    current={}; cof.append(current); continue
                assert current is not None
                cs,es=line.split('|')
                e=tuple(int(z) for z in es.split(','))
                assert len(e)==6 and all(z>=0 for z in e)
                coefficient=M.V.parse_poly(cs,coefficient_ctx,('t',))
                for (power,),v in coefficient.to_dict().items():
                    assert v.q==1 and power<=1
                    key=(power,)+e
                    assert key not in current
                    current[key]=int(v.p)%p
            assert len(cof)==len(polys)
            total=ctx.constant(0)
            for d,g in zip(cof,polys):
                total=reduce(total+ctx.from_dict(d)*g)
            assert total==ctx.constant(1), 'invalid identity '+fn
            print('PASS B0 cell %d roots %d,%d'%(M.CELL,ix,iy),flush=True)
    return True

if __name__=='__main__':
    import cell_model as M
    verify(M)
