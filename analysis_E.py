import json, glob, os, numpy as np
from scipy import stats
R=os.path.join(os.environ.get('V4_ROOT', os.getcwd()), 'results')
def L(e,i): 
    f=f'{R}/{e}/{i}.json'; return json.load(open(f)) if os.path.exists(f) else None
seeds=range(1,11)
def arm(cell, amp, key='near_beta'):
    out=[]
    for s in seeds:
        if cell=='N': r=L('B',f'N_amp{amp}_s{s}')
        elif cell=='P': r=L('B',f'P_amp{amp}_s{s}')
        else: r=L('E',f'{cell}_amp{amp}_s{s}')
        out.append(r[key] if r else np.nan)
    return np.array(out)
def holm(p):
    p=np.asarray(p); o=np.argsort(p); m=len(p); adj=np.empty(m); run=0
    for r,i in enumerate(o):
        run=max(run,(m-r)*p[i]); adj[i]=min(1,run)
    return adj
cells={'N':'normal drive, normal coupling','Ndrive_Pcoup':'normal drive, PD coupling','Pdrive_Ncoup':'PD drive, normal coupling','P':'PD drive, PD coupling'}
print('cell | STN Hz | GPe Hz | beta flat | beta depth1 | gain')
G={}
for c,name in cells.items():
    f=arm(c,0.0); d=arm(c,1.0)
    ok=~np.isnan(f)&~np.isnan(d)
    g=np.exp(np.mean(np.log(d[ok])-np.log(f[ok]))); G[c]=(np.log(d)-np.log(f))
    print(f'{name:32s} {np.nanmean(arm(c,1.0,"stn_rate")):5.1f} {np.nanmean(arm(c,1.0,"gpe_rate")):5.1f}  {np.exp(np.nanmean(np.log(f))):.3g}  {np.exp(np.nanmean(np.log(d))):.3g}  x{g:.2f}  n={ok.sum()}')
# gain contrasts (log gain differences, paired by seed)
rows=[]
for a,b,lab in [('N','Pdrive_Ncoup','drive effect at normal coupling'),('N','Ndrive_Pcoup','coupling effect at normal drive'),
                ('Ndrive_Pcoup','P','drive effect at PD coupling'),('Pdrive_Ncoup','P','coupling effect at PD drive'),('N','P','whole state')]:
    d=G[b]-G[a]; d=d[~np.isnan(d)]
    p=stats.wilcoxon(d).pvalue if len(d)>5 else np.nan
    rows.append((lab,np.exp(d.mean()),p,int((d>0).sum()),len(d)))
ph=holm([r[2] for r in rows])
for r,h in zip(rows,ph): print(f'{r[0]:36s} gain ratio x{r[1]:.2f}  p={r[2]:.3f} Holm={h:.3f}  {r[3]}/{r[4]}')
# gain vs rate across all cells/seeds
rate=np.concatenate([arm(c,1.0,'stn_rate') for c in cells]); g=np.concatenate([G[c] for c in cells])
m=~np.isnan(rate)&~np.isnan(g); print('Spearman gain vs STN rate',stats.spearmanr(rate[m],g[m]))
