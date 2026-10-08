"""Exact necessary mass-ratio restriction on the a=4c branch."""
from pathlib import Path
import json,sympy as s
ROOT=Path(__file__).resolve().parents[1];t,h,d,e,B=s.symbols('t h d e B');a,b,c=s.symbols('a b c')
F=16*t**3-(1+t)**3*(1+t**3);Q=(4*(1+t**3)+1)*(t-1)*(t*t+1)
assert s.monic(s.gcd(F,Q))==t-1
assert s.monic(s.gcd(h**2-4,h**3+8))==h+2
# On this branch the remaining reduced initials imply d^3=1,
# e^2=2*(1+d^2), B=-2e^3, where d=r14/r13,e=r34/r13,B=b/c.
first=d**3-1;second=e**2-2*(1+d*d);third=B+2*e**3
res1=s.factor(s.resultant(first,second,d));assert s.expand(res1-(e*e-4)*(e**4-2*e*e+4))==0
res2=s.factor(s.resultant(res1,third,e));assert s.expand(res2-(B*B-256)*(B*B+32)**2)==0
# Hence positive B implies B=16 (other factors strictly positive).
audit={'a_over_c':4,'forced_t':1,'forced_h':-2,'ratio_constraints':[str(first),str(second),str(third)],'first_resultant':str(res1),'mass_ratio_resultant':str(res2),'positive_mass_necessary_b_over_c':16,'scope':'Necessary conditions only. Full original-input branch derivation audit and remaining equation solvability still required.'}
(ROOT/'results/exceptional_branch_346.json').write_text(json.dumps(audit,indent=2)+'\n');print(audit)
