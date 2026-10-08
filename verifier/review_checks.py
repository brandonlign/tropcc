"""Check elimination input provenance and specialization counterexamples exactly."""
from pathlib import Path
import json
import sympy as sp
import flint
import verify as V

ROOT = Path(__file__).resolve().parents[1]

def main():
    recorded = json.loads((ROOT / 'data/elimination_forms.json').read_text())
    for base, cells in [(334, (334,1767)), (49, (49,2171,2174,2175))]:
        spec = recorded[str(base)]
        for cell in cells:
            header = [f'cell {cell}']
            if cell != 1767:
                header.append('system F+G+H')
            header += ['divmono', 'subset '+','.join(map(str,spec['subset'])),
                       'vars '+','.join(spec['variables'][2:]), '']
            result = V.verify(None, lines=header, ret_gens=True)
            assert len(result) == 4 and result[2] == spec['variables']
            assert list(map(str,result[3][:6])) == spec['polynomials']
        print('PASS elimination inputs', base)
    a,b,x,y = sp.symbols('a b x y')
    for text in recorded['334']['polynomials']:
        poly = sp.Poly(sp.sympify(text.replace('^','**')),a,b)
        assert all(sum(e)==1 for e in poly.monoms())
    print('PASS homogeneous mass reduction (a,b) -> (a/b,1)')
    numerator = (a*y+1)**3
    denominator = a*y**3+1
    assert sp.cancel(denominator.subs(y,-1/a)-(1-1/a**2)) == 0
    assert sp.Poly(numerator,y,domain=sp.QQ.frac_field(a)).gcd(
        sp.Poly(denominator,y,domain=sp.QQ.frac_field(a))).degree() == 0
    phi = numerator-b*b*denominator
    coeffs = sp.Poly(phi,y).all_coeffs()
    assert sp.gcd_list(coeffs) == 1
    print('PASS odd-order nonsquare and primitive phi: irreducibility by Gauss')
    first,second = x+y, x+(1+a)*y+1
    assert sp.expand((second-first).subs(a,0)) == 1
    assert sp.cancel(first.subs({x:1/a,y:-1/a})) == 0
    assert sp.cancel(second.subs({x:1/a,y:-1/a})) == 0
    def initial(expr):
        poly=sp.Poly(expr,x,y)
        degree=max(sum(e) for e in poly.monoms())
        return sum(c*x**e[0]*y**e[1] for e,c in poly.terms() if sum(e)==degree)
    assert initial(first).subs(a,0) == initial(second).subs(a,0) == x+y
    print('PASS negative control: empty special fiber but generic solutions; star detects escape')
    assert all((t*t+1)%3 for t in (1,2))
    q=flint.nmod_poly([1,0,1],3)
    factors=q.factor()[1]
    assert len(factors)==1 and factors[0][0].degree()==2
    t=flint.nmod_poly([0,1],3)
    assert (t*t+1)%q == 0
    print('PASS finite-field negative control: no F_3 torus root, but roots over its closure')
    print('ALL PASS')

if __name__ == '__main__':
    main()
