"""Check every stored maximal cone lies in all 35 tropical hypersurfaces."""
from pathlib import Path
import verify as V
from fan_completeness import maximal_cones

ROOT=Path(__file__).resolve().parents[1]

def main():
    rays,_=V.read_fan(ROOT/'data/fan.out')
    _,polynomials=V.load_sys()
    cones=maximal_cones(ROOT/'data/fan.out')
    for i,f in enumerate(polynomials):
        support=list(f)
        faces=[]
        for ray in rays:
            values=[sum(a*b for a,b in zip(ray,m)) for m in support]
            maximum=max(values)
            faces.append({j for j,v in enumerate(values) if v==maximum})
        for cone in cones:
            common=set.intersection(*(faces[r] for r in cone))
            assert len(common)>=2, (i,cone)
        print('PASS support',i+1,flush=True)
    print('PASS soundness:',len(cones),'maximal cones in all',len(polynomials),'hypersurfaces')

if __name__=='__main__':
    main()
