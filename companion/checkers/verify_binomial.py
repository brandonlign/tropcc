"""Independent direct Laurent-identity verification of exported certificates."""
import hashlib,json,argparse
parser=argparse.ArgumentParser();parser.add_argument('--cells',action='store_true');parser.add_argument('--input');args=parser.parse_args()
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
source=ROOT/(args.input or ('results/cell_binomial_certificates.json' if args.cells else 'results/binomial_certificates.json'))
report=json.loads(source.read_text())
for name,digest in report['input_sha256'].items():
    assert hashlib.sha256((ROOT/'inputs'/name).read_bytes()).hexdigest()==digest
D=json.loads((ROOT/'inputs/sys.json').read_text()); R=s.symbols(D['rvars'])
P=[s.Poly(s.sympify(f),*R) for f in D['polys']]
def monomial(exps): return s.prod(v**e for v,e in zip(R,exps))
def verify(c):
    if not (len(c['generators'])==len(c['shifts'])==len(c['multipliers'])): return False
    if len(c['weight'])!=len(R) or len(c['target_exponent'])!=len(R): return False
    if any(len(shift)!=len(R) for shift in c['shifts']): return False
    if any(not isinstance(i,int) or not 0<=i<len(P) for i in c['generators']): return False
    w=c['weight']; expr=0
    for idx,shift,mult in zip(c['generators'],c['shifts'],c['multipliers']):
        poly=P[idx]; weights={m:sum(e*q for e,q in zip(m,w)) for m in poly.monoms()}
        top=max(weights.values())
        ini=sum(coef*monomial(m) for m,coef in poly.terms() if weights[m]==top)
        expr+=s.sympify(mult)*ini/monomial(shift)
    target=s.sympify(c['parameter_factor'])*monomial(c['target_exponent'])
    return target!=0 and s.expand(expr-target)==0
assert all(verify(c) for c in report['certificates'])
corruption_rejected=None
if report['certificates']:
 corrupt=dict(report['certificates'][0]);corrupt['parameter_factor']='0'
 corruption_rejected=not verify(corrupt)
 assert corruption_rejected
 malformed=dict(report['certificates'][0]);malformed['multipliers']=malformed['multipliers']+['1']
 assert not verify(malformed), 'Mismatched proof arrays must not be silently truncated'
factors=sorted(set(c['parameter_factor'] for c in report['certificates']))
positive_nonzero=all((lambda cs: all(c>0 for c in cs) or all(c<0 for c in cs))(s.Poly(s.sympify(f),*s.symbols('a b c')).coeffs()) for f in factors)
unresolved_factors=[f for f in factors if not (lambda cs: all(c>0 for c in cs) or all(c<0 for c in cs))(s.Poly(s.sympify(f),*s.symbols('a b c')).coeffs())]
out={'certificate_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'verified':len(report['certificates']),'all_factors_nonzero_for_positive_masses':positive_nonzero,'corrupted_certificate_rejected':corruption_rejected,
     'distinct_parameter_factors':factors,'factors_without_coefficient_sign_certificate':unresolved_factors,'sympy':s.__version__,
     'limitations':'Checks identities at stated weights, not fan construction or higher-dimensional cells.'}
(source.with_name(source.stem+'_verification.json') if args.input else ROOT/('results/cell_binomial_verification.json' if args.cells else 'results/binomial_verification.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
