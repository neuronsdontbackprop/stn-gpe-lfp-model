import sys, os, json, glob
ROOT = os.environ.get('V4_ROOT', os.getcwd())
import numpy as np
from scipy import stats
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = ROOT + '/results'
def load(exp):
    return {os.path.basename(f)[:-5]: json.load(open(f)) for f in sorted(glob.glob(f'{R}/{exp}/*.json'))}

def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for r, i in enumerate(o):
        run = max(run, (m - r) * p[i]); adj[i] = min(1, run)
    return adj

def compare(a, b):
    """a, b: paired arrays (same seeds). Log-ratio b/a, Wilcoxon and sign test, both two-sided."""
    d = np.log(b) - np.log(a)
    try: pw = stats.wilcoxon(d, alternative='two-sided').pvalue
    except ValueError: pw = 1.0
    ps = stats.binomtest(int((d > 0).sum()), len(d), 0.5).pvalue
    return dict(ratio=float(np.exp(d.mean())), wilcoxon_p=float(pw), sign_p=float(ps), wins=f'{int((d > 0).sum())}/{len(d)}', n=len(d))

out = {}
F = np.arange(1, 101)

# ------------------------------------------------------------------ B
B = load('B')
AMPS = (0.0, 0.5, 1.0)
if B:
    seeds = sorted({int(k.rsplit('_s', 1)[1]) for k in B})
    def arm(n, a, key, ss=None):
        ss = ss or seeds
        return {s: B[f'{n}_amp{a}_s{s}'][key] for s in ss if f'{n}_amp{a}_s{s}' in B}
    def paired(n1, a1, n2, a2, key):
        x, y = arm(n1, a1, key), arm(n2, a2, key)
        ss = sorted(set(x) & set(y))
        return np.array([x[s] for s in ss]), np.array([y[s] for s in ss])
    res = {}
    for e in ('near', 'far'):
        key = f'{e}_beta'; rows = []
        fam = [('cortex depth 0.5 vs flat (normal)', 'N', 0.0, 'N', 0.5), ('cortex depth 1.0 vs 0.5 (normal)', 'N', 0.5, 'N', 1.0),
               ('cortex depth 0.5 vs flat (PD-like)', 'P', 0.0, 'P', 0.5), ('cortex depth 1.0 vs 0.5 (PD-like)', 'P', 0.5, 'P', 1.0),
               ('state PD-like vs normal, flat cortex', 'N', 0.0, 'P', 0.0), ('state PD-like vs normal, depth 0.5', 'N', 0.5, 'P', 0.5),
               ('state PD-like vs normal, depth 1.0', 'N', 1.0, 'P', 1.0)]
        for name, n1, a1, n2, a2 in fam:
            x, y = paired(n1, a1, n2, a2, key)
            if len(x) >= 4:
                r = compare(x, y); r['name'] = name; rows.append(r)
        # interaction: effect of depth 1.0 vs flat, PD-like vs normal
        pN0, pN1 = paired('N', 0.0, 'N', 1.0, key); pP0, pP1 = paired('P', 0.0, 'P', 1.0, key)
        k = min(len(pN0), len(pP0))
        if k >= 4:
            di = (np.log(pP1[:k]) - np.log(pP0[:k])) - (np.log(pN1[:k]) - np.log(pN0[:k]))
            try: pw = stats.wilcoxon(di).pvalue
            except ValueError: pw = 1.0
            rows.append(dict(name='interaction (state x depth 1.0 vs flat)', ratio=float(np.exp(di.mean())), wilcoxon_p=float(pw),
                             sign_p=float(stats.binomtest(int((di > 0).sum()), len(di), 0.5).pvalue), wins=f'{int((di > 0).sum())}/{len(di)}', n=int(k)))
        if rows:
            ph = holm([r['wilcoxon_p'] for r in rows])
            for r, p in zip(rows, ph): r['holm_p'] = float(p)
        res[e] = rows
    out['B'] = res
    def gm(v):
        v = np.array(list(v)); return float(np.exp(np.log(v).mean())) if len(v) else float('nan')
    out['B_summary'] = {f'{n}_amp{a}': dict(
        n=len(arm(n, a, 'stn_rate')),
        stn=float(np.mean(list(arm(n, a, 'stn_rate').values()))), gpe=float(np.mean(list(arm(n, a, 'gpe_rate').values()))),
        stn_cv=float(np.nanmean(list(arm(n, a, 'stn_cv').values()))), gpe_cv=float(np.nanmean(list(arm(n, a, 'gpe_cv').values()))),
        stn_burst=float(np.mean(list(arm(n, a, 'stn_burst').values()))),
        beta_near=gm(arm(n, a, 'near_beta').values()), beta_far=gm(arm(n, a, 'far_beta').values()),
        dom_share=float(np.mean(list(arm(n, a, 'dom_cell_beta_share').values()))),
        cons=float(np.max(list(arm(n, a, 'cons_ratio').values()))))
        for n in ('N', 'P') for a in AMPS if arm(n, a, 'stn_rate')}
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    labels = [(n, a) for n in ('N', 'P') for a in AMPS]
    for axx, e in zip(ax[:2], ('near', 'far')):
        for i, (n, a) in enumerate(labels):
            v = np.array(list(arm(n, a, f'{e}_beta').values()))
            if len(v) == 0: continue
            axx.scatter(np.full(len(v), i) + np.random.default_rng(0).uniform(-.12, .12, len(v)), v, s=16, alpha=.65)
            axx.plot([i - .25, i + .25], [np.exp(np.log(v).mean())] * 2, 'k-', lw=2)
        axx.set_yscale('log'); axx.set_xticks(range(len(labels))); axx.set_xticklabels([f'{"normal" if n=="N" else "PD-like"}\ndepth {a}' for n, a in labels], fontsize=8)
        axx.set_title(f'LFP beta power (13-30 Hz), {e} electrode\n(each dot = one seed, bar = geometric mean)', fontsize=9); axx.set_ylabel('uV^2')
    for i, (n, a) in enumerate(labels):
        s = list(arm(n, a, 'stn_rate').values()); g = list(arm(n, a, 'gpe_rate').values())
        if s: ax[2].bar(i - .2, np.mean(s), .38, color='C0'); ax[2].bar(i + .2, np.mean(g), .38, color='C1')
    ax[2].axhline(13.5, color='C0', ls=':'); ax[2].axhline(30.9, color='C0', ls=':'); ax[2].axhline(32.0, color='C1', ls=':'); ax[2].axhline(17.4, color='C1', ls=':')
    ax[2].set_xticks(range(len(labels))); ax[2].set_xticklabels([f'{"N" if n=="N" else "PD"}\n{a}' for n, a in labels], fontsize=8)
    ax[2].set_title('Mean rates: STN (blue) and GPe (orange), Hz\ndotted = Cruz et al. 2011 control / 6-OHDA', fontsize=9)
    plt.tight_layout(); plt.savefig(ROOT + '/fig_B_states_cortex.png', dpi=130); plt.close()

    # spectra (seed-averaged) for the report
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for n, ls in (('N', '-'), ('P', '--')):
        for a, c in zip(AMPS, ('gray', 'C0', 'C3')):
            v = arm(n, a, 'near_psd')
            if v:
                P = np.mean(list(v.values()), 0); ax[0].semilogy(F, P, ls, color=c, label=f'{"normal" if n=="N" else "PD-like"}, depth {a}')
            v = arm(n, a, 'stn_pop_psd')
            if v:
                P = np.mean(list(v.values()), 0); ax[1].semilogy(F, P, ls, color=c)
    ax[0].axvspan(13, 30, color='y', alpha=.15); ax[1].axvspan(13, 30, color='y', alpha=.15)
    ax[0].set_title('LFP spectrum (near electrode), seed-averaged', fontsize=9); ax[1].set_title('STN population-rate spectrum', fontsize=9)
    ax[0].set_xlabel('Hz'); ax[1].set_xlabel('Hz'); ax[0].legend(fontsize=7)
    plt.tight_layout(); plt.savefig(ROOT + '/fig_B_spectra.png', dpi=130); plt.close()

# ------------------------------------------------------------------ A
A = load('A')
Ds = (3, 5, 10, 15); Ks = (1, 3, 6)
if A:
    def seeds_of(d, k): return [s for s in (1, 2, 3, 4) if f'd{d}_k{k}_s{s}' in A]
    def agg(key, d, k):
        v = [A[f'd{d}_k{k}_s{s}'][key] for s in seeds_of(d, k)]; return np.mean(v) if v else np.nan
    def pk(d, k, lo=5, hi=60, which='stn_pop_psd'):
        ss = seeds_of(d, k)
        if not ss: return np.nan
        P = np.mean([A[f'd{d}_k{k}_s{s}'][which] for s in ss], 0); m = (F >= lo) & (F <= hi)
        return float(F[m][np.argmax(P[m])])
    def prom(d, k):
        ss = seeds_of(d, k)
        if not ss: return np.nan
        P = np.mean([A[f'd{d}_k{k}_s{s}']['stn_pop_psd'] for s in ss], 0); m = (F >= 13) & (F <= 30)
        return float(P[m].max() / np.median(P))
    grid = {key: np.array([[agg(key, d, k) for d in Ds] for k in Ks]) for key in ('stn_rate', 'gpe_rate', 'near_beta', 'near_tot')}
    grid['bfrac'] = grid['near_beta'] / grid['near_tot']
    grid['pop_peak'] = np.array([[pk(d, k) for d in Ds] for k in Ks]); grid['prom'] = np.array([[prom(d, k) for d in Ds] for k in Ks])
    out['A'] = {k: v.tolist() for k, v in grid.items()}; out['A']['n_runs'] = len(A)
    fig, ax = plt.subplots(1, 4, figsize=(17, 3.6))
    for axx, (key, title, cm) in zip(ax, [('stn_rate', 'STN rate (Hz)', 'viridis'), ('bfrac', 'LFP beta fraction (13-30 / 1-100 Hz)', 'magma'),
                                          ('prom', 'STN pop-rate: peak in 13-30 Hz / median', 'magma'), ('pop_peak', 'STN pop-rate peak frequency (Hz, 5-60)', 'cividis')]):
        im = axx.imshow(grid[key], aspect='auto', cmap=cm, origin='lower')
        axx.set_xticks(range(len(Ds))); axx.set_xticklabels(Ds); axx.set_yticks(range(len(Ks))); axx.set_yticklabels([f'x{k}' for k in Ks])
        axx.set_xlabel('delay each way (ms)'); axx.set_ylabel('coupling scale (w_sg, w_gs)'); axx.set_title(title, fontsize=9)
        for i in range(len(Ks)):
            for j in range(len(Ds)):
                axx.text(j, i, f'{grid[key][i, j]:.2g}', ha='center', va='center', color='w', fontsize=7)
        plt.colorbar(im, ax=axx, fraction=0.046)
    plt.tight_layout(); plt.savefig(ROOT + '/fig_A_loop_only.png', dpi=130); plt.close()

# ------------------------------------------------------------------ C
C = load('C')
if C:
    names = ['base', 'gpe_up', 'gpe_down', 'no_STN_to_GPe', 'no_GPe_to_STN', 'no_loop']
    rowsC = {}
    for nm in names:
        ks = [k for k in C if k.rsplit('_s', 1)[0] == nm]
        if not ks: continue
        def pkf(which, lo=35, hi=70):
            P = np.mean([C[k][which] for k in ks], 0); m = (F >= lo) & (F <= hi)
            return float(F[m][np.argmax(P[m])]), float(P[m].max() / np.median(P))
        rowsC[nm] = dict(n=len(ks), stn=float(np.mean([C[k]['stn_rate'] for k in ks])), gpe=float(np.mean([C[k]['gpe_rate'] for k in ks])),
                         stn_pop_peak_35_70=pkf('stn_pop_psd'), lfp_peak_35_70=pkf('near_psd'), gpe_pop_peak_35_70=pkf('gpe_pop_psd'))
    out['C'] = rowsC


# ------------------------------------------------------------------ D  (exploratory: is the ~35 Hz bump of the PD-like state a loop rhythm?)
D = load('D')
if D and B:
    s5 = [s for s in range(1, 6)]
    def pack(label, recs):
        P = np.mean([r['stn_pop_psd'] for r in recs], 0); L = np.mean([r['near_psd'] for r in recs], 0); m = (F >= 22) & (F <= 45)
        b = lambda X, lo, hi: float(X[(F >= lo) & (F <= hi)].sum())
        return dict(label=label, n=len(recs), stn_rate=float(np.mean([r['stn_rate'] for r in recs])), stn_cv=float(np.nanmean([r['stn_cv'] for r in recs])),
                    gpe_rate=float(np.mean([r['gpe_rate'] for r in recs])), peak_22_45=float(F[m][np.argmax(P[m])]),
                    prominence=float(P[m].max() / np.median(P)), stn_pop_30_40=b(P, 30, 40), lfp_30_40=b(L, 30, 40), spectrum=P.tolist())
    rowsD = [pack('normal', [B[f'N_amp0.0_s{s}'] for s in s5 if f'N_amp0.0_s{s}' in B]),
             pack('PD-like (full loop)', [B[f'P_amp0.0_s{s}'] for s in s5 if f'P_amp0.0_s{s}' in B])]
    for v, lab in (('PD_no_STN_to_GPe', 'PD-like, no STN->GPe'), ('PD_no_GPe_to_STN', 'PD-like, no GPe->STN'), ('PD_no_loop', 'PD-like, no loop')):
        rr = [D[f'{v}_s{s}'] for s in s5 if f'{v}_s{s}' in D]
        if rr: rowsD.append(pack(lab, rr))
    out['D'] = rowsD
    # 10-seed paired test of the bump, normal vs PD-like (exploratory: bands chosen after looking at the spectra)
    if all(f'{n}_amp0.0_s{s}' in B for n in 'NP' for s in range(1, 11)):
        ex = {}
        for which in ('stn_pop_psd', 'near_psd'):
            for lo, hi in ((13, 30), (30, 40), (40, 60)):
                m = (F >= lo) & (F <= hi)
                x = np.array([np.sum(np.array(B[f'N_amp0.0_s{s}'][which])[m]) for s in range(1, 11)])
                y = np.array([np.sum(np.array(B[f'P_amp0.0_s{s}'][which])[m]) for s in range(1, 11)])
                ex[f'{which}_{lo}-{hi}Hz'] = compare(x, y)
        ph = holm([v['wilcoxon_p'] for v in ex.values()])
        for v, p in zip(ex.values(), ph): v['holm_p'] = float(p)
        out['B_flat_bands_exploratory'] = ex
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.8))
    for r, c in zip(rowsD, ('gray', 'C3', 'C1', 'C2', 'C0')):
        ax.semilogy(F, r['spectrum'], color=c, label=f"{r['label']} (STN {r['stn_rate']:.0f} Hz)")
    ax.axvspan(13, 30, color='y', alpha=.15); ax.set_xlabel('Hz'); ax.set_title('STN population-rate spectrum, flat cortex, 5 seeds\nthe 30-40 Hz bump survives cutting the loop', fontsize=9)
    ax.legend(fontsize=7); plt.tight_layout(); plt.savefig(ROOT + '/fig_D_bump.png', dpi=130); plt.close()

json.dump(out, open(ROOT + '/summary4.json', 'w'), indent=1)
print('saved summary4.json', {k: (len(v) if hasattr(v, '__len__') else v) for k, v in out.items()})
