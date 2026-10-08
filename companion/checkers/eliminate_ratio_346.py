"""Necessary low-dimensional conditions from cell 346, with exact checks."""
from pathlib import Path
import json,sympy as s
ROOT=Path(__file__).resolve().parents[1];a,b,c=s.symbols('a b c');t,h,z=s.symbols('t h z')
# t=r25/r15, h=r12/r15, z=r15^3. All are nonzero.
# H from subsystem and CM generator 32 give:
H=c*h**3+a*(1+t**3)
C=h**2*t-(1+t)*(1+t**3)
# Generator 26, after substitution and removal of invertible factors:
Q=(a*(1+t**3)+c)*(h**2*t-2-t-t**3)
Qmod=s.rem(Q,s.Poly(C,h),h).as_expr()
assert s.expand(Qmod-(a*(1+t**3)+c)*(t**4-1))==0
F=s.factor(s.resultant(H,C,h))
expected=a**2*t**3-c**2*(1+t)**3*(1+t**3)
# Resultant has an extra factor (1+t^3)^2; nonzero by H in the torus.
assert s.expand(F-(1+t**3)**2*expected)==0
B=(a*(1+t**3)+c)*(t-1)*(t**2+1)
# Q mod C = (t+1)*B. t+1 nonzero follows from 1+t^3 != 0.
assert s.expand(Qmod-(t+1)*B)==0
res=s.factor(s.resultant(expected,B,t))
assert res!=0
out={'cell':346,'ratio_constraints':[str(H),str(C),str(Q)],'eliminated_ratio_equations':[str(expected),str(B)],'necessary_mass_polynomial':str(res),'sample_nonzero_at_a_1_c_1':str(res.subs({a:1,c:1})),'scope':'Necessary polynomial obstruction derived from selected generators via torus divisions. Independent full certificate and branch audit still required before adding any exclusions.'}
(ROOT/'results/ratio_346.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
