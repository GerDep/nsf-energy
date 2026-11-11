"""Re-express archived NSF integration samples in the attached variance note's notation.

No amplitude or integration sample is changed. E_i=1 throughout.
The code's old Q is the note's W; C is accumulated as (E-1)^2 P.
"""
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
WORK=ROOT.parent/'reference_data'
INPUTS=set()

def read(rel):
    p=WORK/rel;INPUTS.add(p)
    if p.suffix=='.jsonl':return [json.loads(l) for l in p.read_text().splitlines()]
    return json.loads(p.read_text())

def basis(r):
    return [r['P'],(r['E']-1)*r['P'],(r['E']-1)**2*r['P'],r['E']**2*r['P']]

def reduce23(rows,Nlow=256,Nuv=128,R=8,high=None,upper=None):
    values=[];metadata=[]
    for band,N in [('low',Nlow),('uv',Nuv)]:
        seeds=sorted({r['seed'] for r in rows if r['band']==band});assert len(seeds)==R
        vv=[]
        for seed in seeds:
            rr=sorted((r for r in rows if r['band']==band and r['seed']==seed),key=lambda r:r['index'])
            assert len(rr)==N and [r['index'] for r in rr]==list(range(N))
            vectors=[basis(r) if upper is None or r['E']<=upper else [0.]*4 for r in rr]
            if high is not None:
                hh=[r for r in high if r['band']==band and r['seed']==seed]
                lut={r['index']:r for r in rr}
                assert len(hh)==len({r['index'] for r in hh})
                for r in hh:
                    b=lut[r['index']];assert r['E']==b['E'] and r['P_low']==b['P']
                    x,y=basis(r),basis(b);p=r['inclusion_probability'];assert 0<p<=1
                    vectors.append([(a-c)/p for a,c in zip(x,y)])
            v=[math.fsum(x[j] for x in vectors)/N for j in range(4)]
            vv.append(v);metadata.append(dict(band=band,seed=seed,N=N,estimate_PDCQ=v,
                                             estimate_PBCQ=[v[0],v[0]+v[1],v[2],v[3]]))
        values.append(np.array(vv))
    return values,metadata

def fold22(rows,maximum=2024.):
    values=[]
    for r in rows:
        u,P,D=r['u'],r['P'],r['D']
        if D>=0:C=(u-1)*D;Q=u*u*P
        elif u<=maximum:C=-2*D;Q=P
        else:C=-(u-1)/u*D;Q=P/(u*u)
        values.append([P,D,C,Q])
    return np.array([math.fsum(x[j] for x in values) for j in range(4)])

def variance(v):
    P,D,C=v[:3];return C/P-(D/P)**2

def grad_v(v):
    P,D,C=v[:3];return np.array([-C/P**2+2*D**2/P**3,-2*D/P**2,1/P])

def covariance(groups):
    return sum(np.cov(g[:,:3],rowvar=False,ddof=1)/len(g) for g in groups)

def error(groups,g):
    g=np.asarray(g,dtype=np.longdouble);v=np.longdouble(0)
    for a in groups:
        x=np.asarray(a[:,:len(g)],dtype=np.longdouble);d=(x-x.mean(0))@g
        v+=np.sum(d*d)/(len(x)*(len(x)-1))
    return float(np.sqrt(v))

def summarize(m2,groups):
    m3=sum(g.mean(0) for g in groups);total=m2+m3;P,D,C,Q=total;r=D/P;v=variance(total)
    assert P>0 and v>0
    cov=covariance(groups);L=np.array([[1,0,0],[1,1,0],[0,0,1]])
    cov_PBC=L@cov@L.T
    gv=grad_v(total);g_note=np.array([-C/P**2+2*r*(P+D)/P**2,-2*r/P,1/P])
    ev=error(groups,gv);ev_note=math.sqrt(g_note@cov_PBC@g_note)
    assert math.isclose(ev,ev_note,rel_tol=1e-9)
    channels=[]
    for n,m in [(2,m2),(3,m3)]:
        p,d,c,q=m;rn=d/p;vn=variance(m);assert p>0 and vn>0
        channels.append(dict(channel=n,P=float(p),W=float(p+d),Q=float(q),C=float(c),
                             weight=float(p/P),mu=float(1+rn),r_E=float(rn),v=float(vn),
                             mu_SE=0. if n==2 else error(groups,[-d/p**2,1/p,0]),
                             v_SE=0. if n==2 else error(groups,grad_v(m)),
                             deterministic=(n==2),Q_identity_residual=float(q-(p+2*d+c))))
    a,b=channels;within=a['weight']*a['v']+b['weight']*b['v']
    between=a['weight']*b['weight']*(a['r_E']-b['r_E'])**2
    rawv=Q/P-((P+D)/P)**2
    checks=dict(C_identity_absolute=float(abs(C-(Q-2*(P+D)+P))),
                raw_vs_centered_variance=float(abs(rawv-v)),
                mixture_variance_absolute=float(abs(within+between-v)),
                positivity_margin=float(C/P-r*r),note_gradient_relative_SE_difference=ev_note/ev-1)
    assert checks['raw_vs_centered_variance']<1e-12 and checks['mixture_variance_absolute']<1e-12
    assert all(abs(c['Q_identity_residual'])<1e-10*max(1.,c['P']) for c in channels)
    jackvar=0.
    for a in groups:
        xs=np.array([variance(total+(a.mean(0)-row)/(len(a)-1)) for row in a])
        jackvar+=(len(a)-1)/len(a)*np.sum((xs-xs.mean())**2)
    return dict(channels=channels,P=float(P),W=float(P+D),Q=float(Q),C=float(C),r_E=float(r),mu=float(1+r),
                v=float(v),v_SE=ev,sd=float(np.sqrt(v)),sd_SE=float(ev/(2*np.sqrt(v))),
                mean_SE=error(groups,[-D/P**2,1/P,0]),mean_estimator_variance=error(groups,[-D/P**2,1/P,0])**2,
                weight_SE=error(groups,[-m2[0]/P**2,0,0]),within=float(within),between=float(between),
                covariance_PDC=cov.tolist(),covariance_PBC=cov_PBC.tolist(),
                v_jackknife_SE=float(np.sqrt(jackvar)),checks=checks)

def direct_cut(tag):
    rows=read(f'{tag}/muestras23_n16.jsonl');summary=read(f'{tag}/summary.json')
    groups,meta=reduce23(rows)
    fine=fold22(read(f'{tag}/cuadratura22_e6_q6_n16.jsonl'))
    coarse=fold22(read(f'{tag}/cuadratura22_e4_q6_n16.jsonl'))
    s=summarize(fine,groups);s.update(kappa=summary['kappa'],Emax=2024.,R=8,Nlow=256,Nuv=128,
        Ntotal=3072,profile='direct refined angular order 16',replicas=meta,
        radial_v_change=float(variance(fine+sum(g.mean(0) for g in groups))-variance(coarse+sum(g.mean(0) for g in groups))),
        angular_control_nodes=summary['counts']['angular_checks'],angular_control_max=summary['angular_control_max_relative_change'])
    assert math.isclose(s['r_E'],(summary['D22']+summary['D23'])/(summary['I22']+summary['I23']),rel_tol=1e-12)
    assert math.isclose(s['mean_SE'],summary['combined_mean_MC_SE'],rel_tol=1e-10)
    return s,groups

def old_cut():
    z=read('figure_sources/calculo_IR_1e-5.json');raw=z['raw_samples'];s=z['summary']
    groups,meta=reduce23(raw['deep_base.jsonl'],R=4,high=raw['deep_high_n32_fine.jsonl'])
    fine=fold22(raw['deep22_k1e-05_n16_e6_q6_seg1_fine.jsonl'])
    coarse=fold22(raw['deep22_k1e-05_n16_e4_q6_seg1_fine.jsonl'])
    v=sum(g.mean(0) for g in groups)
    out=summarize(fine,groups);out.update(kappa=1e-5,Emax=2024.,R=4,Nlow=256,Nuv=128,Ntotal=1536,
       profile='order-8 base plus inclusion-weighted refined order-32 correction',replicas=meta,
       refined_nodes=len(raw['deep_high_n32_fine.jsonl']),radial_v_change=float(variance(fine+v)-variance(coarse+v)),
       angular_control_max=s['angular_validation_max_relative_change_23'])
    assert math.isclose(out['mu'],s['combined_mean'],abs_tol=1e-14)
    assert math.isclose(out['mean_SE'],s['combined_mean_MC_SE'],rel_tol=1e-10)
    return out

def upper_scan(base_groups):
    tail=read('paper_nsf_audit/uv_tail_n16.jsonl');old=read('paper_nsf_audit/emax_complete.json')
    assert old['complete'] and old['validation_passed']
    results=[];allgroups=[];jac=[]
    for j,M in enumerate([2024.,4048.,16192.,18096.]):
        # Tail has a single band and eight independent seeds.
        vals=[]
        for seed in sorted({r['seed'] for r in tail}):
            rr=sorted((r for r in tail if r['seed']==seed),key=lambda r:r['index']);assert len(rr)==128
            vals.append([math.fsum(basis(r)[k] for r in rr if r['E']<=M)/128 for k in range(4)])
        groups=base_groups+[np.array(vals)]
        if M==2024.:
            fine=fold22(read('calculo_IR_1e-8/cuadratura22_e6_q6_n16.jsonl'),M)
            coarse=fold22(read('calculo_IR_1e-8/cuadratura22_e4_q6_n16.jsonl'),M)
        else:
            fine=fold22(read(f'paper_nsf_audit/quad22_E{int(M)}_ne6.jsonl'),M)
            coarse=fold22(read(f'paper_nsf_audit/quad22_E{int(M)}_ne4.jsonl'),M)
        s=summarize(fine,groups);s.update(kappa=1e-8,Emax=M,Ntotal=3072+(1024 if M>2024 else 0),
            radial_v_change=float(variance(fine+sum(g.mean(0) for g in groups))-variance(coarse+sum(g.mean(0) for g in groups))))
        assert math.isclose(s['r_E'],old['points'][j]['r_E'],rel_tol=1e-12)
        g=np.zeros(12);g[j*3:j*3+3]=grad_v(fine+sum(x.mean(0) for x in groups));jac.append(g)
        allgroups.append(groups);results.append(s)
    joint=[np.concatenate([allgroups[j][b][:,:3] for j in range(4)],axis=1) for b in range(3)]
    for j,s in enumerate(results):
        s['delta_v_vs_2024']=s['v']-results[0]['v'];s['delta_v_SE']=error(joint,jac[j]-jac[0])
    return dict(points=results,angular_control_nodes=len(old['angular']),angular_control_max=old['max_angular_relative_change'])

def angular_partition():
    config=read('angular_2048/comparacion_radio_angular_IR1e-7_Emax2048.json');assert config['validation_passed']
    base=read('angular_2048/muestras23_n16.jsonl');half=read('angular_2048/muestras23_mitad_n16.jsonl')
    lut={r['id']:r for r in half};groups=[]
    for band,N in [('low',256),('uv',128)]:
        vals=[]
        for seed in sorted({r['seed'] for r in base if r['band']==band}):
            rr=[r for r in base if r['band']==band and r['seed']==seed];assert len(rr)==N
            vectors=[]
            for r in rr:
                h=lut[r['id']];assert r['E']==h['E'] and r['u5']==h['u5']
                dp=h['P']-r['P'];d=r['E']-1
                vectors.append([r['P'],d*r['P'],d*d*r['P'],dp,d*dp,d*d*dp])
            vals.append([math.fsum(v[j] for v in vectors)/N for j in range(6)])
        groups.append(np.asarray(vals,dtype=np.longdouble))
    m=sum(g.mean(0) for g in groups)
    bf=fold22(read('angular_2048/cuadratura22_e6_q6_n16.jsonl'),2048.)[:3]
    hf=fold22(read('angular_2048/cuadratura22_mitad_e6_q6_n16.jsonl'),2048.)[:3]
    b=m[:3]+bf;h=b+m[3:]+(hf-bf);P,D,C=b;dp,dd,dc=h-b
    rb=D/P;rh=h[1]/h[0]
    delta=(dc*P-C*dp)/(P*(P+dp))-(rh-rb)*(rh+rb)
    gb=np.r_[grad_v(b),np.zeros(3)];gh=np.r_[grad_v(h),grad_v(h)]
    return dict(kappa=1e-7,Emax=2048.,base_fraction=.35,half_fraction=.175,base_v=float(variance(b)),
       half_v=float(variance(h)),delta_v=float(delta),paired_SE=error(groups,gh-gb),
       base_v_SE=error(groups,gb),angular_control_nodes=82,angular_control_max=config['angular_controls']['max_absolute_relative_P_change'])

def main():
    global WORK
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-dir',type=Path,default=ROOT.parent/'reference_data')
    ap.add_argument('--output',type=Path,default=ROOT.parent/'results/generated/variance.json')
    args=ap.parse_args()
    WORK=args.data_dir.expanduser().resolve()
    required=WORK/'calculo_IR_1e-8/muestras23_n16.jsonl'
    if not required.is_file():
        ap.error('Reference samples not found. First run: python scripts/import_reference_data.py nsf_reference_data.zip')
    INPUTS.clear()
    small=[];basegroups=None
    for tag in ['calculo_IR_1e-6','calculo_IR_1e-7','calculo_IR_1e-8']:
        s,g=direct_cut(tag);small.append(s)
        if tag.endswith('1e-8'):basegroups=g
    rows=[old_cut()]+small
    upper=upper_scan(basegroups);angular=angular_partition()
    result=dict(convention='Ei=1; note P, W, Q=second raw moment, C=centered about Ei; old code Q means W',
        complete=True,rows=rows,upper_energy=upper,angular_partition=angular,
        source_note='Calculo_varianza_energia_NSF.pdf, 11 November 2026',
        inputs={str(p.relative_to(WORK)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(INPUTS)})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for s in rows:
        print(s['kappa'],'p2,p3',[c['weight'] for c in s['channels']],'mu2,mu3',[c['mu'] for c in s['channels']],
              'v2,v3',[c['v'] for c in s['channels']],'V +/-',s['v'],s['v_SE'],'radial',s['radial_v_change'])
    print('UPPER:',[(s['Emax'],s['v'],s['v_SE'],s['delta_v_vs_2024'],s['delta_v_SE']) for s in upper['points']])
    print('ANGULAR:',angular)

if __name__=='__main__':main()
