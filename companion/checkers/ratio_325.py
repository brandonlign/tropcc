"""Exact necessary positive-mass ratio at cell 325 from three initials."""
from pathlib import Path
import json,sympy as s,hashlib
ROOT=Path(__file__).resolve().parents[1];D=json.loads((ROOT/'inputs/sys.json').read_text());R=s.symbols(D['rvars']);a,b,c=s.symbols('a b c');text=(ROOT/'inputs/fan.out').read_text()
def block(k):return text.split('\n'+k+'\n',1)[1].split('\n\n',1)[0].strip().splitlines()
rays=[tuple(map(int,row.split('#')[0].split())) for row in block('RAYS')];cones=[tuple(map(int,row.split('#')[0].strip().strip('{}').split())) for row in block('CONES')];cones=[C for C in cones if C];w=tuple(sum(rays[i][j] for i in cones[325]) for j in range(10))
def initial(i):
 terms=s.Poly(s.sympify(D['polys'][i]),*R).terms();scores=[sum(x*y for x,y in zip(m,w)) for m,coef in terms];top=max(scores);terms=[row for row,score in zip(terms,scores) if score==top];shift=tuple(min(m[j] for m,coef in terms) for j in range(10));return s.expand(sum(coef*s.prod(x**(e-d) for x,e,d in zip(R,m,shift)) for m,coef in terms))
T,U,A=s.symbols('T U A');x=R[1];sub={a:A*b,R[2]:T*x,R[7]:U*x}
f=1+T*T-U*U+2*T**3;g=2+T+T**3-T*U*U;h=-2*A*U**3-1+T*T-U*U
for i,expected,scale in [(1,f,-b*x**3),(2,g,-b*x**3),(8,h,b*x**3)]:assert s.cancel(initial(i).subs(sub,simultaneous=True)/scale-expected)==0
assert s.expand(g-T*f-2*(1-T**4))==0
# T=0 is impossible; split the exact factors of T^4-1.
assert s.factor(T**4-1)==(T-1)*(T+1)*(T*T+1)
# T=-1 forces U^2=0, excluded in torus.
assert s.expand(f.subs(T,-1)+U*U)==0
# T=1: U^2=4, h=0 gives AU^3=-2, hence 16 A^2=1.
p=s.Poly(h.subs(T,1),U);h1=p.rem(s.Poly(U*U-4,U)).as_expr();assert s.expand(h1+8*A*U+4)==0
assert s.expand(s.resultant(U*U-4,2*A*U+1,U)-(1-16*A*A))==0
# T^2=-1: f implies U^2=-2T; h reduces to 4AUT+2(T-1).
G=s.groebner([T*T+1,U*U+2*T],U,T,domain=s.QQ.frac_field(A));assert s.expand(G.reduce(h)[1]-(4*A*U*T+2*T-2))==0
# Squaring 2AUT=1-T and using those relations gives 4 A^2+1=0.
assert G.reduce((2*A*U*T)**2-(1-T)**2-2*T*(4*A*A+1))[1]==0
out={'cell':325,'used_generators':[1,2,8],'ratio_variables':{'T':'r14/r13','U':'r34/r13','A':'a/b'},'necessary_mass_polynomial':'(16*a**2-b**2)*(4*a**2+b**2)','positive_mass_necessary_condition':'b=4*a','positive_mass_forced_distance_ratios':{'r14/r13':1,'r34/r13':-2},'scope':'Necessary conditions only; exceptional b=4a not resolved by this subsystem. No full torus solution claimed.','input_sha256':{f:hashlib.sha256((ROOT/'inputs'/f).read_bytes()).hexdigest() for f in ['sys.json','fan.out']}}
(ROOT/'results/ratio_325.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
