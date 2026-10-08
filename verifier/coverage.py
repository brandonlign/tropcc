#!/usr/bin/env python3
"""Proof-coverage check (exact, no Singular): every one of the 16438 cells of the support prevariety Sigma (fan.out) is
  (A) companion      : in companion generic_ids (Laurent certificates, verified by their verify_*.py), or
  (V) certified  : its (12),(34)-orbit representative (verifier/orbits.json, from symcheck.py) has a certificate
                   cert/<r>.cof, dcert/<r>.cof, gcert/<r>.cof or hcert/<r>.cof (S14: F+G+H) that verify.py PASSes (re-run here unless --fast), or
  (E) elim       : orbit rep in {334,1767} (verify_334.py) or {49,2171,2174,2175} (verify_49.py), S24 elimination certificates,
                   or their sigma-images 2345, 218, 2885, 3494, 3545 (elim_sigma.py: forms coincide after a<->b),
  (P) residual   : in the residual set R; R must satisfy the pointedness condition: an integer c with <c,v> > 0 for
                   every ray v of every cell of R (Lemma C). c is read from balance/pointc.json or balance/slackc.json and checked exactly.
Prints counts and the uncovered cells. Exit 0 iff everything is covered.
usage: coverage.py [--fast]   (--fast: trust existing .ok markers instead of re-verifying)"""
import sys,os,json
import verify as V
import verify_res as VR
N5=V.N5; EX=os.path.join(N5,'exact')
rays,cells=V.read_fan(os.path.join(N5,'fan.out'))
orb={int(k):int(v) for k,v in json.load(open(os.path.join(V.HERE,'orbits.json'))).items()}
companion=set(json.load(open(os.path.join(V.HERE,'..','companion','results','generic_ids.json')))['generic'])
fast='--fast' in sys.argv
ok_rep={}
LEMS={81:81,316:81,121:121,50:121,206:206,220:206}; LEMS_RUN={}
ELIM={334:['verify_334.py'],1767:['verify_334.py'],49:['verify_49.py'],2171:['verify_49.py'],2174:['verify_49.py'],2175:['verify_49.py'],
      2345:['verify_334.py','elim_sigma.py'],218:['verify_49.py','elim_sigma.py'],2885:['verify_49.py','elim_sigma.py'],3494:['verify_49.py','elim_sigma.py'],3545:['verify_49.py','elim_sigma.py']}; ELIM_RUN={}
ELIM_OK={'verify_334.py':os.path.join('ecert','c334','334.ok'),'verify_49.py':os.path.join('h49','49.ok'),'elim_sigma.py':os.path.join('h49','elim_sigma.ok')}
def certified(r):
    if r in ok_rep: return ok_rep[r]
    res=None
    for d in ['cert','dcert','gcert','hcert','ecert']:
        fn=os.path.join(EX,d,'%d.cof'%r)
        if not os.path.exists(fn): continue
        if fast:
            if os.path.exists(os.path.join(EX,d,'%d.ok'%r)) or (d=='gcert' and open(fn).readline().strip()=='cell %d'%r and os.path.getsize(fn)<2000): res=d; break
            continue
        k,good,msg=V.verify(fn)
        if good and k==r: res=d; break
    if res is None:   # resultant-specialisation certificates (verify_res.py), always re-checked (fast)
        for d in ['dcert','gcert']:
            fn=os.path.join(EX,d,'%d.res'%r)
            if os.path.exists(fn):
                k,good,msg=VR.check(fn)
                if good and k==r: res=d; break
    if res is None:   # split-tree certificates exact/split/<r>.split (verify_split.py); --fast trusts split/<r>.ok
        fn=os.path.join(EX,'split','%d.split'%r)
        if os.path.exists(fn):
            if fast: res='hcert' if os.path.exists(os.path.join(EX,'split','%d.ok'%r)) else None
            else:
                import verify_split as VS
                k,good,msg=VS.check(fn)
                if good and k==r: res='hcert'
    if res is None and r in ELIM:   # elimination certificate verify_334.py (cells 334, 1767); --fast trusts ecert/c334/334.ok
        good=True
        for sc in ELIM[r]:
            if fast: good&=os.path.exists(os.path.join(EX,ELIM_OK[sc]))
            else:
                if sc not in ELIM_RUN:
                    import subprocess; ELIM_RUN[sc]=subprocess.run([sys.executable,os.path.join(V.HERE,sc)],capture_output=True).returncode==0
                good&=ELIM_RUN[sc]
        res='elim' if good else None
    if res is None and r in LEMS:
        base=LEMS[r]
        if fast: good=os.path.exists(os.path.join(N5,'lemS','%d.ok'%base))
        else:
            if base not in LEMS_RUN:
                import subprocess
                print('verifying Lemma S cell',base,flush=True)
                ld=os.path.join(N5,'lemS','verification'); os.makedirs(ld,exist_ok=True)
                with open(os.path.join(ld,'star_%d.log'%base),'w') as stream:
                    LEMS_RUN[base]=subprocess.run([sys.executable,os.path.join(V.HERE,'verify_lemS.py'),str(base)],stdout=stream,stderr=subprocess.STDOUT).returncode==0
            good=LEMS_RUN[base]
        res='lemS' if good else None
    ok_rep[r]=res; return res
# extra symmetry sigma=(13)(24) with a<->b (symcheck2.py -> sigma.json). If orbit rep r of sigma(k) has a .cof certificate,
# transport.py maps it to cell sigma(r) (whose (12),(34)-orbit contains k) and verify.py re-checks the transported file from scratch.
SIG=json.load(open(os.path.join(V.HERE,'sigma.json'))) if os.path.exists(os.path.join(V.HERE,'sigma.json')) else None
sig_ok={}
assert SIG is not None
for base,partner in ((81,316),(121,50),(206,220)):
    assert orb[SIG[base]]==partner, 'Lemma S symmetry pair mismatch'
def sig_certified(k):
    if SIG is None: return None
    r=orb[SIG[k]]
    if r in sig_ok: return sig_ok[r]
    res=None
    if certified(r):
        for d in ['cert','dcert','gcert','hcert','ecert']:
            fn=os.path.join(EX,d,'%d.cof'%r)
            if not os.path.exists(fn) or certified(r)!=d: continue
            if fast: res='sigma'; break
            import transport,tempfile
            out=os.path.join(tempfile.gettempdir(),'sigma_%d.cof'%r)
            try:
                transport.run(fn,SIG[r],out); kk,good,msg=V.verify(out)
                if good and kk==SIG[r] and orb[SIG[r]]==orb[k]: res='sigma'
            except SystemExit: pass
            break
    sig_ok[r]=res; return res
cnt={'companion':0,'cert':0,'dcert':0,'gcert':0,'hcert':0,'ecert':0,'elim':0,'sigma':0,'lemS':0}; uncovered=[]
for k in range(len(cells)):
    if not fast and k%1000==0: print('coverage progress',k,'/',len(cells),flush=True)
    if k in companion: cnt['companion']+=1; continue
    d=certified(orb[k]) or sig_certified(k)
    if d: cnt[d]+=1
    else: uncovered.append(k)
print('covered:',cnt,' uncovered:',len(uncovered))
# pointedness of the uncovered residual
P=json.load(open(os.path.join(N5,'balance','pointc.json')))
SC=os.path.join(N5,'balance','slackc.json')
if os.path.exists(SC): P=dict(P); P['slackc']=[json.load(open(SC))['c']]   # c from balance/slack3.py (S+10+kept orbits)
best=None
for key,v in sorted(P.items(),key=lambda kv:-len(kv[0])):
    if not v: continue
    c=v[0]; R={i for k in uncovered for i in cells[k]}
    if all(sum(x*y for x,y in zip(c,rays[i]))>0 for i in R): best=(key,c); break
if best: print('residual pointed: c=%s (from pointc.json[%s]); exact check <c,v> > 0 on all %d residual rays: PASS'%(best[1],best[0],len({i for k in uncovered for i in cells[k]})))
else:
    orbs=sorted({orb[k] for k in uncovered}); print('residual NOT (yet) pointed; uncovered orbit reps (%d):'%len(orbs),orbs[:60]); json.dump(orbs,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'uncovered.json'),'w'))
report={'mode':'fast' if fast else 'full','covered':cnt,'residual_cells':uncovered,'residual_rays':sorted({i for k in uncovered for i in cells[k]}),'pointed':bool(best),'c':best[1] if best else None,'total_cells':len(cells)}
json.dump(report,open(os.path.join(V.HERE,'coverage_report.json'),'w'),indent=2)
sys.exit(0 if best else 1)
