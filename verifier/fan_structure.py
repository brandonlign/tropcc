"""Verify the stored fan as cones of the common Newton normal fan.

Exact rational facet enumeration, without gfan, cddlib, or SMT. Each cone
must equal its exposed-face normal region and every facet must be stored.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
import argparse, hashlib, json, math
import flint
import verify as V
from fan_completeness import check_symmetries

ROOT = Path(__file__).resolve().parents[1]

def primitive(row, canonical=False):
    denominator=math.lcm(*(int(x.denominator) for x in row))
    integers=[int(x*denominator) for x in row]
    divisor=math.gcd(*integers)
    if not divisor:
        return None
    integers=[x//divisor for x in integers]
    if canonical and next(x for x in integers if x)<0:
        integers=[-x for x in integers]
    return tuple(integers)

def geometry(rays, cone, supports, faces):
    equalities=set(); inequalities=set()
    for support,rayfaces in zip(supports,faces):
        face=set.intersection(*(rayfaces[r] for r in cone))
        assert len(face)>=2, ('not tropical',cone)
        base=support[min(face)]
        for j,m in enumerate(support):
            row=primitive([a-b for a,b in zip(base,m)],j in face)
            if row:
                (equalities if j in face else inequalities).add(row)
    basis=V.nullspace_Q(sorted(equalities),10)
    d=len(basis); assert 1<=d<=4
    # Select an invertible coordinate submatrix of the basis.
    matrix=flint.fmpq_mat(d,10,[x for b in basis for x in b])
    reduced,rank=matrix.rref(); assert rank==d
    pivots=[next(j for j in range(10) if reduced[i,j]) for i in range(d)]
    inverse=flint.fmpq_mat(d,d,[basis[j][i] for i in pivots for j in range(d)]).inv()
    coords=[]
    for r in cone:
        solution=inverse*flint.fmpq_mat(d,1,[rays[r][i] for i in pivots])
        t=[solution[i,0] for i in range(d)]
        assert all(sum(basis[j][i]*t[j] for j in range(d))==rays[r][i] for i in range(10))
        coords.append(primitive(t))
    assert flint.fmpq_mat(len(coords),d,[x for t in coords for x in t]).rank()==d
    normals=set()
    for selected in combinations(coords,d-1):
        null=V.nullspace_Q(selected,d)
        if len(null)!=1:
            continue
        h=primitive(null[0]); values=[V.dot(h,t) for t in coords]
        if all(v<=0 for v in values):
            h=tuple(-x for x in h); values=[-v for v in values]
        if all(v>=0 for v in values) and any(values):
            normals.add(h)
    assert normals
    restricted={primitive([V.dot(row,b) for b in basis]):row for row in sorted(inequalities)}
    # C is contained in N by its maximizing faces; every facet of C is a
    # defining inequality of N, proving the reverse containment directly.
    assert normals<=restricted.keys(), ('missing defining facet',cone,normals-restricted.keys())
    facets=[]
    for h in sorted(normals):
        facets.append(tuple(sorted(r for r,t in zip(cone,coords) if V.dot(h,t)==0)))
    for r,t in zip(cone,coords):
        assert sum(V.dot(h,t) for h in normals)>0, ('lineality',r)
        active=[h for h in normals if V.dot(h,t)==0]
        assert [s for s,u in zip(cone,coords) if all(V.dot(h,u)==0 for h in active)]==[r], ('nonextreme ray',r)
    halfspaces=[(True,row) for row in sorted(equalities)]
    halfspaces += [(False,restricted[h]) for h in sorted(normals)]
    return d,facets,halfspaces

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit',type=int,help='Development subset; not a full audit')
    args=parser.parse_args()
    rays,cells=V.read_fan(ROOT/'data/fan.out')
    cone_block=(ROOT/'data/fan.out').read_text().split('\nCONES\n',1)[1].split('\n\n',1)[0]
    assert any(line.split('#')[0].strip()=='{}' for line in cone_block.splitlines()), 'missing origin'
    rv,polynomials=V.load_sys()
    assert len(set(map(tuple,rays)))==len(rays)
    assert all(primitive(r)==tuple(r) for r in rays)
    stored={tuple(sorted(c)) for c in cells}; assert len(stored)==len(cells)
    check_symmetries(rays,cells,rv,polynomials)
    supports=[list(f) for f in polynomials]; faces=[]
    for support in supports:
        rayfaces=[]
        for ray in rays:
            values=[V.dot(ray,m) for m in support]; maximum=max(values)
            rayfaces.append({j for j,v in enumerate(values) if v==maximum})
        faces.append(rayfaces)
    dimensions=Counter(); facets=0; descriptions=[]
    for k,cone in enumerate(cells[:args.limit]):
        d,boundary,h=geometry(rays,cone,supports,faces)
        assert all(not f or f in stored for f in boundary), ('missing face',k,boundary)
        dimensions[d]+=1; facets+=len(boundary); descriptions.append(h)
        if k%1000==0:
            print('Checked',k,'cones',flush=True)
    report={'complete':args.limit is None,'nonzero_cones':sum(dimensions.values()),
            'dimensions':dict(sorted(dimensions.items())),'facets_checked':facets,
            'normal_regions_equal':True,'closed_under_faces':True,
            'rays_extreme_and_primitive':True,'full_fan_symmetry_checked':True,
            'fan_sha256':hashlib.sha256((ROOT/'data/fan.out').read_bytes()).hexdigest(),
            'system_sha256':hashlib.sha256((ROOT/'data/sys.json').read_bytes()).hexdigest()}
    if args.limit is None:
        (ROOT/'data/fan_structure_audit.json').write_text(json.dumps(report,indent=2)+'\n')
        (ROOT/'logs/fan_exact_halfspaces.json').write_text(json.dumps(descriptions)+'\n')
    print('PASS',report,flush=True)

if __name__=='__main__':
    main()
