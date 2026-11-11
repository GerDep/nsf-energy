"""Create the variance-note tables from the checked moment reduction."""
import argparse, json, math
from pathlib import Path

HERE=Path(__file__).resolve().parent
WORK=HERE.parent
OUT=WORK/'results/generated'

def sci(x,d=3):
    if not x:return '0'
    e=math.floor(math.log10(abs(x)));a=x/10**e
    return rf'{a:.{d}f}\times10^{{{e}}}'

def pm(x,se,d=3):
    if not se:return sci(x,d)
    e=math.floor(math.log10(abs(x)));scale=10.**e
    return rf'({x/scale:.{d}f}\pm{se/scale:.{d}f})\times10^{{{e}}}'

def fixed_pm(x,se):
    decimals=max(0,1-math.floor(math.log10(se)))
    return rf'{x:.{decimals}f}\pm{se:.{decimals}f}'

def mathcell(x):return '$'+x+'$'
def row(cells):return ' & '.join(cells)+r'\\'

def main():
    global OUT
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=WORK/'results/reference/variance.json')
    ap.add_argument('--output-dir',type=Path,default=WORK/'results/generated')
    args=ap.parse_args()
    OUT=args.output_dir.expanduser().resolve()
    r=json.loads(args.input.read_text());assert r['complete']
    rows=r['rows'];markers={};channel=[];variances=[];moments=[];parts=[];sampling=[];radial=[]
    for s in rows:
        k=mathcell(rf'10^{{{round(math.log10(s["kappa"]))}}}')
        c2,c3=s['channels'];pe=100*s['weight_SE']
        channel.append(row([k,mathcell(fixed_pm(100*c2['weight'],pe)),mathcell(fixed_pm(100*c3['weight'],pe)),
                            mathcell(f'{c2["mu"]:.12f}'),mathcell(f'{c3["mu"]:.12f}')]))
        variances.append(row([k,mathcell(sci(c2['v'])),mathcell(sci(c3['v'])),
                             mathcell(pm(s['v'],s['v_SE'])),mathcell(pm(s['sd'],s['sd_SE']))]))
        for c in s['channels']:
            moments.append(row([k,str(c['channel'])]+[mathcell(sci(c[n],7)) for n in ['P','W','Q','C']]))
        parts.append(row([k]+[mathcell(sci(x,7)) for x in [s['within'],s['between'],s['r_E']**2]]))
        sampling.append(row([k,str(s['R']),str(s['Ntotal']),mathcell(sci(s['r_E'],5)),mathcell(sci(s['mean_SE'],4))]))
        radial.append(row([k,mathcell(sci(abs(s['radial_v_change'])))]))
    markers.update(CHANNEL_ROWS='\n'.join(channel),VARIANCE_ROWS='\n'.join(variances),MOMENT_ROWS='\n'.join(moments),
                   DECOMPOSITION_ROWS='\n'.join(parts),SAMPLING_ROWS='\n'.join(sampling),RADIAL_ROWS='\n'.join(radial))
    upper=r['upper_energy']['points'];markers['UPPER_ROWS']='\n'.join(row([f'{s["Emax"]:.0f}',mathcell(pm(s['v'],s['v_SE'])),
        mathcell(pm(s['sd'],s['sd_SE'])),mathcell(sci(s['r_E'],4))]) for s in upper)
    ang=r['angular_partition']
    markers.update(RAW_ERROR=sci(max(s['checks']['raw_vs_centered_variance'] for s in rows)),
       MIXTURE_ERROR=sci(max(s['checks']['mixture_variance_absolute'] for s in rows)),
       JACKKNIFE_MAX=f"{max(abs(s['v_jackknife_SE']/s['v_SE']-1) for s in rows)*100:.2f}",
       UV_RATIO=f"{upper[-1]['v']/upper[0]['v']:.2f}",UV_DELTA=pm(upper[-1]['delta_v_vs_2024'],upper[-1]['delta_v_SE']),
       ANGLE_BASE=sci(ang['base_v'],9),ANGLE_HALF=sci(ang['half_v'],9),ANGLE_DELTA=sci(ang['delta_v']),
       ANGLE_DELTA_SE=sci(ang['paired_SE']),ANGLE_BASE_SE=sci(ang['base_v_SE']))
    tex=(WORK/'docs/variance_report_template.tex').read_text()
    for k,v in markers.items():tex=tex.replace('@@'+k+'@@',v)
    assert '@@' not in tex
    OUT.mkdir(exist_ok=True,parents=True)
    path=OUT/'varianza_NSF_formato_nota.tex';path.write_text(tex)
    lines=['DATOS NSF EN EL FORMATO DE LA NOTA DE VARIANZA','Ei=1; Emax/Ei=2024; errores estadisticos de una sigma.',
           'P=normalizacion; W=primer momento; Q=segundo momento; C=segundo momento centrado en Ei.',
           'Todos los resultados corresponden a las amplitudes retenidas.','',
           'kappa/Ei | n | pn | mu_n | var_n | Pn | Wn | Qn | Cn']
    for s in rows:
        for c in s['channels']:
            lines.append(' | '.join([f'{s["kappa"]:g}',str(c['channel'])]+[f'{c[k]:.16g}' for k in ['weight','mu','v','P','W','Q','C']]))
    lines+=['','kappa/Ei | var_total | error_var | sd_fisica | error_sd | error_media | replicas | puntos']
    for s in rows:lines.append(' | '.join([f'{s["kappa"]:g}']+[f'{s[k]:.12e}' for k in ['v','v_SE','sd','sd_SE','mean_SE']]+[str(s['R']),str(s['Ntotal'])]))
    lines+=['','SENSIBILIDAD A EMAX, kappa/Ei=1e-8','Emax/Ei | var_total | error_var | sd_fisica | delta_var_vs_2024 | error_delta_pareado']
    for s in upper:lines.append(' | '.join([f'{s["Emax"]:g}']+[f'{s[k]:.12e}' for k in ['v','v_SE','sd','delta_v_vs_2024','delta_v_SE']]))
    (OUT/'tablas_varianza_NSF.txt').write_text('\n'.join(lines)+'\n')
    print(path)
    print(OUT/'tablas_varianza_NSF.txt')

if __name__=='__main__':main()
