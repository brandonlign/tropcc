"""Independent original-input audit of the ratio obstruction at cell 346."""
from pathlib import Path
import json,sympy as s
ROOT=Path(__file__).resolve().parents[1];D=json.loads((ROOT/'inputs/sys.json').read_text());R=s.symbols(D['rvars']);a,b,c=s.symbols('a b c');p,x,y,v,u,k,q,z,l,n=R
text=(ROOT/'inputs/fan.out').read_text()
def block(key):return text.split('\n'+key+'\n',1)[1].split('\n\n',1)[0].strip().splitlines()
rays=[tuple(map(int,row.split('#')[0].split())) for row in block('RAYS')];cones=[tuple(map(int,row.split('#')[0].strip().strip('{}').split())) for row in block('CONES')];cones=[C for C in cones if C];w=tuple(sum(rays[i][j] for i in cones[346]) for j in range(10))
g=[]
for f in D['polys']:
 poly=s.Poly(s.sympify(f),*R);terms=poly.terms();top=max(sum(i*j for i,j in zip(m,w)) for m,coef in terms);terms=[(m,coef) for m,coef in terms if sum(i*j for i,j in zip(m,w))==top];shift=[min(m[j] for m,coef in terms) for j in range(10)];g.append(s.expand(sum(coef*s.prod(X**(e-d) for X,e,d in zip(R,m,shift)) for m,coef in terms)))
H0=a*(v**3+q**3)+c*p**3
J1=v**3*u**2-(v**3+q**3)*l**2
J2=v**3*k**2-(v**3+q**3)*n**2
assert s.expand(g[1]-g[18]-H0*l**2)==0
assert s.expand(g[18]-a*J1)==0
assert s.expand(g[19]-a*J2)==0
# Reduce by proved consequences over the field of remaining variables.
K=s.QQ.frac_field(a,b,c,x,y,v,q,z,l,n)
G=s.groebner([H0,J1,J2],p,u,k,domain=K)
def rem(f):return G.reduce(s.expand(f))[1]
t,h=s.symbols('t h');H=a*(1+t**3)+c*h**3;C=h*h*t-(1+t)*(1+t**3);Q=(a*(1+t**3)+c)*(h*h*t-2-t-t**3)
sub={t:q/v,h:p/v}
assert s.cancel(rem(g[32])-2*q**2*n**4*C.subs(sub))==0
assert s.cancel(rem(g[26])-(-a*v**6/c)*Q.subs(sub))==0
assert s.cancel(H.subs(sub)-H0/v**3)==0
F=a*a*t**3-c*c*(1+t)**3*(1+t**3);B=(a*(1+t**3)+c)*(t-1)*(t*t+1)
assert s.expand(Q-(a*(1+t**3)+c)*C-(t+1)*B)==0
assert s.expand(s.resultant(H,C,h)-(1+t**3)**2*F)==0
# Verify nonzero resultant also by direct Sylvester determinant at a=c=1.
# This suffices to check generic nonvanishing independently, not every coefficient.
f=s.Poly(F.subs({a:1,c:1}),t);bb=s.Poly(B.subs({a:1,c:1}),t)
nf,nb=f.degree(),bb.degree();rows=[]
for coeffs,number in [(f.all_coeffs(),nb),(bb.all_coeffs(),nf)]:
 for shift in range(number):rows.append([0]*shift+coeffs+[0]*(nf+nb-shift-len(coeffs)))
determinant=s.Matrix(rows).det();assert determinant==50625
mass=s.sympify(json.loads((ROOT/'results/ratio_346.json').read_text())['necessary_mass_polynomial']);assert s.expand(s.resultant(F,B,t)-mass)==0
out={'cell':346,'original_input_initials_recomputed':True,'subsystem_consequences_verified':True,'CM_and_AC_remainders_verified':True,'ratio_resultant_verified':True,'independent_specialized_sylvester_determinant':str(determinant),'necessary_mass_polynomial':str(mass),'nonvanishing_assumptions':['a','c','all ten distances'],'derived_nonvanishing':['1+t^3 (from H and nonzero c,h)','1+t (divides 1+t^3)'],'conclusion':'The original initial-generator system at cell 346 has no torus solution over Q(a,b,c); for fixed masses the nonzero resultant is a sufficient exclusion condition with a*c != 0.'}
(ROOT/'results/ratio_346_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
