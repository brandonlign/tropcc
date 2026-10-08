# Independent exact polynomial I/O (no sympy, no Singular). Polynomials = dict {exponent tuple: int coeff}.
import re
def parse_expanded(s, names):
    """Parse an expanded polynomial string (sympy str or Singular string) into {exp tuple: int}.
    Accepts '**' or '^' powers, integer coefficients, '*' products; names = ordered variable list."""
    idx={v:i for i,v in enumerate(names)}; n=len(names)
    s=s.replace('**','^').replace(' ','')
    if s[0] not in '+-': s='+'+s
    out={}
    for sign,body in re.findall(r'([+-])([^+-]+)',s):
        c=1; e=[0]*n
        for f in body.split('*'):
            if f=='' : raise ValueError(s)
            if f[0].isdigit():
                if '/' in f: raise ValueError('rational coefficient not allowed: '+f)
                c*=int(f); continue
            if '^' in f: v,p=f.split('^'); p=int(p)
            else: v,p=f,1
            e[idx[v]]+=p
        c=-c if sign=='-' else c
        t=tuple(e); out[t]=out.get(t,0)+c
        if out[t]==0: del out[t]
    return out
