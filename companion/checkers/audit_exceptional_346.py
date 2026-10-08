"""Check branch reductions against fresh original initial equations."""
import runpy,json
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
# Recompute and audit originals; obtain independently reconstructed generators.
ctx=runpy.run_path(str(ROOT/'checkers/audit_ratio_346.py'));g=ctx['g'];R=ctx['R'];a,b,c=ctx['a'],ctx['b'],ctx['c'];p,x,y,v,u,k,q,z,l,n=R
branch={a:4*c,p:-2*v,q:v};gg=[s.expand(f.subs(branch,simultaneous=True)) for f in g]
# Previously proved subsystem ratios give u^2=2l^2,k^2=2n^2.
K=s.QQ.frac_field(b,c,x,y,v,z,l,n)
Gb=s.groebner([u*u-2*l*l,k*k-2*n*n],u,k,domain=K)
r=[s.expand(Gb.reduce(f)[1]) for f in gg]
assert s.expand(r[16]-b*v**4*(l*l+n*n))==0
# b,v nonzero, so n^2=-l^2.
r=[s.rem(f,s.Poly(n*n+l*l,n),n).as_expr() for f in r]
assert s.expand(r[3]-b*l*l*(y**3-x**3))==0
assert s.expand(r[11]+2*l*l*(b*x**3+2*c*z**3))==0
assert s.expand(r[31]+2*l**4*(2*x*x+2*y*y-z*z))==0
# Positive b/c and the exact resultant imply b=16c.
# z^3=-8x^3. Setting e=z/x, d=y/x gives d^3=1,e^2=2(1+d^2).
d,e=s.symbols('d e');GB=s.groebner([d**3-1,e*e-2*(1+d*d),e**3+8],e,d,domain=s.QQ)
assert GB.reduce(e+2)[1]==0 and GB.reduce(d-1)[1]==0
# Thus y=x,z=-2x. Generator 8 is a nonzero torus monomial.
last=s.factor(r[8].subs({b:16*c,y:x,z:-2*x},simultaneous=True))
assert s.expand(last+72*c*x**4*l**2)==0
out={'cell':346,'original_initial_branch_reductions_verified':True,'positive_mass_exceptional_ratios':{'a/c':4,'b/c':16},'forced_distance_ratios':{'r14/r13':1,'r34/r13':-2},'final_generator_8':str(last),'conclusion':'Combined with audited ratio resultant and Sturm argument, cell 346 is excluded for ALL positive a,b,c. No torus solution of its initial-generator system exists for those masses.','scope':'One cell only; uses exact elimination and branch reasoning, not yet a single exported Laurent identity.'}
(ROOT/'results/exceptional_346_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
