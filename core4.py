# -*- coding: utf-8 -*-
"""
v4 core: STN-GPe network whose STN cells are the PUBLISHED rat subthalamic
projection neuron of Gillies & Willshaw (2006, J Neurophysiol 95:2352-2365;
ModelDB 74298), run at 35 degC.  The LFP is computed from the membrane
currents of exactly these simulated cells (line-source approximation).

What is NOT published here: the 3-D layout of the dendrites (the original model
has only lengths/diameters, so branch directions are drawn at random for each
cell), the GPe cell (a reduced Hodgkin-Huxley point cell tuned in v3, with its
rates referenced to 35 degC so it keeps its v3 behaviour), the synapses and the
network wiring.  All of that is stated in the report.
"""
import os, platform, subprocess
import numpy as np
from scipy.signal import welch
from scipy.spatial.transform import Rotation as Rot
from neuron import h
h.load_file('stdrun.hoc')

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------- loading
def load_all():
    """Compile/load channel files and the published cell template (once)."""
    mod_dir = os.path.join(HERE, 'mods')
    arch = 'x86_64' if platform.machine() == 'x86_64' else 'aarch64'
    lib = os.path.join(mod_dir, arch, 'libnrnmech.so')
    if not os.path.exists(lib):
        subprocess.run(['nrnivmodl'], cwd=mod_dir, capture_output=True, text=True, check=True)
    try:
        h.nrn_load_dll(lib)
    except RuntimeError:
        pass                                  # already loaded (e.g. auto-loaded)
    cwd = os.getcwd(); os.chdir(HERE)          # SThprotocell.hoc reads sth-data/ relatively
    try:
        h.xopen(os.path.join(HERE, 'SThprotocell.hoc').replace(os.sep, '/'))
        proto = h.SThproto()
    finally:
        os.chdir(cwd)
    return proto

_PROTO = None
def proto():
    global _PROTO
    if _PROTO is None:
        _PROTO = load_all()
    return _PROTO

# conductance maps of the published cell (soma + dendrites), read once
GMAPS = {}
def gmap(name):
    if name not in GMAPS:
        GMAPS[name] = np.loadtxt(os.path.join(HERE, 'sth-data', f'cell-{name}'))
    return GMAPS[name]

CONDUCTANCES = ['gk_KDR', 'gk_Kv31', 'gk_Ih', 'gk_sKCa', 'gcaT_CaT', 'gcaN_HVA', 'gcaL_HVA']

# -------------------------------------------------------------------- config
CFG = dict(
    T_SIM=3000.0, T_DISCARD=1000.0, DT=0.05, REC_DT=0.25,
    CELSIUS=35.0,
    N_STN=20, N_GPE=20,
    POP_RADIUS=150.0,
    ELEC_NEAR=(0., 0., 200.), ELEC_FAR=(0., 0., 800.),
    SIGMA=0.3e-3,           # 0.3 S/m in nA/(uV*um)
    P_CONN=0.4,
)

DEFAULT = dict(
    w_gs=0.0012, w_sg=0.0012, w_gg=0.0006,       # uS
    d_sg=5.0, d_gs=5.0, d_gg=2.0,                # ms
    stn_bias=0.0, stn_bias_sd=0.0,               # nA, mean and cell-to-cell sd
    bg_rate=0.0, bg_w=0.0,                       # background (non-cortical) AMPA Poisson: Hz, uS
    ctx_mode='poisson', ctx_base=150.0, ctx_w=0.0170, ctx_amp=0.0, ctx_freq=20.0,
    ctx_delay=4.0, ctx_phase_jitter=False,
    gpe_bias=0.0, gpe_bias_sd=0.0,               # nA into each GPe cell (mean, cell-to-cell sd)
)

# --------------------------------------------------------------- build cells
def build_stn(idx):
    """Create one published STN cell and set the paper's conductances."""
    cell = h.SThcell(idx, proto())
    soma = cell.soma
    n0, n1 = int(proto().ntree0), int(proto().ntree1)
    secs = [soma] + [cell.dend0[i] for i in range(n0)] + [cell.dend1[i] for i in range(n1)]
    soma.gna_Na = 1.483419823e-02
    soma.gna_NaL = 1.108670852e-05
    for nm in CONDUCTANCES:
        for tree, sec, x, val in gmap(nm):
            tree = int(tree); sec = int(sec)
            seg = soma(0.5) if tree == -1 else (cell.dend0[sec](x) if tree == 0 else cell.dend1[sec](x))
            setattr(seg, nm, float(val))
    for s in secs:                       # ion styles exactly as in sample.hoc
        h.ion_style('na_ion', 1, 2, 1, 0, 1, sec=s)
        h.ion_style('k_ion', 1, 2, 1, 0, 1, sec=s)
        h.ion_style('ca_ion', 3, 2, 1, 1, 1, sec=s)
    return cell, secs

def layout_3d(secs, rng):
    """Random branch directions, children start where NEURON attaches them
    (parent_x along the parent), so LFP geometry == simulated topology."""
    soma = secs[0]; L0 = soma.L; diam0 = soma.diam
    geom = {soma.name(): (np.array([0., 0., -L0 / 2]), np.array([0., 0., L0 / 2]))}
    dirs = {soma.name(): np.array([0., 0., 1.])}
    queue = [soma]
    while queue:
        par = queue.pop(0)
        for ch in par.children():
            px = ch.parentseg().x
            ps, pe = geom[par.name()]
            a = ps + px * (pe - ps)
            if par.name() == soma.name():
                base = np.array([0., 0., 1.0 if px > 0.5 else -1.0])
            else:
                base = dirs[par.name()]
            v = base + 0.9 * rng.normal(0, 1, 3)
            v /= np.linalg.norm(v)
            geom[ch.name()] = (a, a + ch.L * v); dirs[ch.name()] = v
            queue.append(ch)
    for s in secs:
        a, b = geom[s.name()]
        d = s.diam
        s.pt3dclear(); s.pt3dadd(*a, d); s.pt3dadd(*b, d)
    return geom

def build_gpe(prefix):
    s = h.Section(name=f'{prefix}_soma')
    s.L = 15.0; s.diam = 15.0; s.Ra = 123.0; s.cm = 1.0; s.nseg = 1
    s.insert('gpehh'); s.gnabar_gpehh = 0.20; s.gkbar_gpehh = 0.022
    s.gl_gpehh = 0.0007; s.el_gpehh = -62.0
    s.insert('gpecat'); s.ggpecatbar_gpecat = 0.0005; s.eca = 120.0
    return s

# ----------------------------------------------------------- LFP forward model
def seg_table(secs, geom, cell_idx):
    rows = []
    for s in secs:
        a, b = geom[s.name()]
        n = s.nseg
        for i in range(n):
            rows.append((cell_idx, a + (i / n) * (b - a), a + ((i + 1) / n) * (b - a), s.diam / 2))
    return rows

def line_kernel(ws, we, radius, elec, sigma):
    ds = we - ws; L = float(np.linalg.norm(ds))
    if L < 1e-9:
        d = max(float(np.linalg.norm(0.5 * (ws + we) - elec)), radius)
        return 1.0 / (4 * np.pi * sigma * d)
    u = ds / L; d = elec - ws
    hp = float(d @ u); rp = max(float(np.linalg.norm(d - hp * u)), radius)
    return (np.arcsinh(hp / rp) - np.arcsinh((hp - L) / rp)) / (4 * np.pi * sigma * L)

def positions(seed):
    rng = np.random.default_rng(10_000 + seed)
    n = CFG['N_STN']
    u = rng.uniform(0, 1, n) ** (1 / 3)
    th = np.arccos(rng.uniform(-1, 1, n)); ph = rng.uniform(0, 2 * np.pi, n)
    R0 = CFG['POP_RADIUS']
    pos = np.stack([R0 * u * np.sin(th) * np.cos(ph), R0 * u * np.sin(th) * np.sin(ph), R0 * u * np.cos(th)], 1)
    ori = Rot.random(n, random_state=10_000 + seed).as_matrix()
    return pos, ori

# ------------------------------------------------------------------ cortex
def cortical_train(seed, t_stop, base_hz, beta_hz, beta_amp, phase=0.0):
    rng = np.random.default_rng(seed)
    tg = np.arange(0, t_stop + 1.0, 1.0)
    rate = np.clip(base_hz * (1 + beta_amp * np.cos(2 * np.pi * beta_hz * tg / 1000.0 + phase)), 0, None)
    rmax = max(base_hz * (1 + beta_amp) * 1.05, 1e-9)
    n = rng.poisson(rmax * t_stop / 1000.0)
    cand = np.sort(rng.uniform(0, t_stop, n))
    keep = rng.uniform(0, 1, n) < np.interp(cand, tg, rate) / rmax
    return cand[keep]

def poisson_train(seed, t_stop, rate_hz):
    rng = np.random.default_rng(seed)
    n = rng.poisson(rate_hz * t_stop / 1000.0)
    return np.sort(rng.uniform(0, t_stop, n))

# ---------------------------------------------------------------- analysis
def spikes_from_v(v, t, thr=-10.0, t0=0.0):
    above = v > thr
    idx = np.where((~above[:-1]) & above[1:])[0]
    st = t[idx]
    return st[st > t0]

def spike_stats(spk, T, t0):
    rate = len(spk) / ((T - t0) / 1000.0)
    if len(spk) > 3:
        isi = np.diff(spk); isi = isi[isi > 0]
        return rate, float(isi.std(ddof=1) / isi.mean()), float(np.mean(isi < 10.0))
    return rate, float('nan'), 0.0

def band(f, P, lo, hi):
    m = (f >= lo) & (f <= hi)
    return float(np.trapezoid(P[m], f[m]))

BANDS = dict(dt=(1, 8), alpha=(8, 13), beta=(13, 30), lgam=(30, 60), hi=(60, 100), tot=(1, 100))

def psd(x, fs):
    f, P = welch(x - np.mean(x), fs=fs, nperseg=int(fs))
    return f, P

# -------------------------------------------------------------- simulation
def simulate(seed, lfp=True, **over):
    """One network run. Call in a fresh process (published cells are hoc objects)."""
    p = dict(DEFAULT); p.update(over)
    proto()
    h.cvode.use_fast_imem(1)
    rng = np.random.default_rng(30_000 + seed)
    nS, nG = CFG['N_STN'], CFG['N_GPE']
    pos, ori = positions(seed)

    cells, stn_secs, geoms = [], [], []
    for i in range(nS):
        c, secs = build_stn(i)
        g = layout_3d(secs, np.random.default_rng(20_000 + 100 * seed + i))
        cells.append(c); stn_secs.append(secs); geoms.append(g)
    gpe = [build_gpe(f'gpe{i}') for i in range(nG)]
    # aCSF of Bevan & Wilson (1999), as in the published sample file
    h.nai0_na_ion = 15; h.nao0_na_ion = 128.5; h.ki0_k_ion = 140; h.ko0_k_ion = 2.5
    h.cai0_ca_ion = 1e-4; h.cao0_ca_ion = 2.0

    # synapses: GABA on a proximal dendrite, AMPA (cortex/background) on a distal terminal
    syn_ctx, syn_gaba_s, ics, ics_stn = [], [], [], []
    for i in range(nS):
        c = cells[i]
        a = h.ExpSyn(c.dend0[7](0.5)); a.tau = 2.0; a.e = 0.0; syn_ctx.append(a)
        g = h.ExpSyn(c.dend0[1](0.5)); g.tau = 6.0; g.e = -80.0; syn_gaba_s.append(g)
        amp = p['stn_bias'] + p['stn_bias_sd'] * rng.normal()
        if amp != 0.0:
            ic = h.IClamp(c.soma(0.5)); ic.delay = 0; ic.dur = 1e9; ic.amp = amp; ics.append(ic); ics_stn.append(ic)
    syn_ampa_g, syn_gaba_g = [], []
    for j in range(nG):
        e = h.ExpSyn(gpe[j](0.5)); e.tau = 2.0; e.e = 0.0; syn_ampa_g.append(e)
        g = h.ExpSyn(gpe[j](0.5)); g.tau = 6.0; g.e = -80.0; syn_gaba_g.append(g)
        gb = p['gpe_bias'] + p['gpe_bias_sd'] * rng.normal()
        if gb != 0.0:
            ic = h.IClamp(gpe[j](0.5)); ic.delay = 0; ic.dur = 1e9; ic.amp = gb; ics.append(ic)

    Wsg = rng.uniform(0, 1, (nG, nS)) < CFG['P_CONN']
    Wgs = rng.uniform(0, 1, (nS, nG)) < CFG['P_CONN']
    Wgg = rng.uniform(0, 1, (nG, nG)) < CFG['P_CONN']; np.fill_diagonal(Wgg, False)
    ncs = []
    def wire(src_sec, syn, w, d):
        nc = h.NetCon(src_sec(0.5)._ref_v, syn, sec=src_sec)
        nc.threshold = -10.0; nc.delay = d; nc.weight[0] = w; ncs.append(nc)
    for j in range(nG):
        for i in range(nS):
            if Wsg[j, i]: wire(cells[i].soma, syn_ampa_g[j], p['w_sg'], p['d_sg'])
    for i in range(nS):
        for j in range(nG):
            if Wgs[i, j]: wire(gpe[j], syn_gaba_s[i], p['w_gs'], p['d_gs'])
    for j in range(nG):
        for k in range(nG):
            if Wgg[j, k]: wire(gpe[k], syn_gaba_g[j], p['w_gg'], p['d_gg'])

    queues = []
    for i in range(nS):
        if p['ctx_mode'] == 'poisson':
            ph = float(rng.uniform(0, 2 * np.pi)) if p['ctx_phase_jitter'] else 0.0
            sp = cortical_train(50_000 + 1000 * seed + i, CFG['T_SIM'], p['ctx_base'], p['ctx_freq'], p['ctx_amp'], ph)
            nc = h.NetCon(None, syn_ctx[i]); nc.weight[0] = p['ctx_w']
            queues.append((nc, sp + p['ctx_delay']))
        if p['bg_rate'] > 0:
            sp = poisson_train(70_000 + 1000 * seed + i, CFG['T_SIM'], p['bg_rate'])
            nc = h.NetCon(None, syn_ctx[i]); nc.weight[0] = p['bg_w']
            queues.append((nc, sp))
    def _load():
        for nc, ts in queues:
            for t in ts: nc.event(float(t))
    fih = h.FInitializeHandler(_load)

    vS = [h.Vector().record(cells[i].soma(0.5)._ref_v, 0.1) for i in range(nS)]
    vG = [h.Vector().record(gpe[j](0.5)._ref_v, 0.1) for j in range(nG)]
    tv = h.Vector().record(h._ref_t, 0.1)
    imem, rows, cell_slices = [], [], []
    if lfp:
        for i in range(nS):
            start = len(imem)
            for s in stn_secs[i]:
                for seg in s:
                    imem.append(h.Vector().record(seg._ref_i_membrane_, CFG['REC_DT']))
            rows += seg_table(stn_secs[i], geoms[i], i)
            cell_slices.append((start, len(imem)))
        ti = h.Vector().record(h._ref_t, CFG['REC_DT'])

    h.dt = CFG['DT']; h.celsius = CFG['CELSIUS']
    h.v_init = -65.0
    h.finitialize(-65.0); h.continuerun(CFG['T_SIM'])

    T0 = CFG['T_DISCARD']
    t = np.array(tv)
    out = dict(seed=seed, params={k: (v if not isinstance(v, np.generic) else v.item()) for k, v in p.items()})
    spS = [spikes_from_v(np.array(v), t, t0=T0) for v in vS]
    spG = [spikes_from_v(np.array(v), t, t0=T0) for v in vG]
    stS = [spike_stats(s, CFG['T_SIM'], T0) for s in spS]; stG = [spike_stats(s, CFG['T_SIM'], T0) for s in spG]
    out['stn_rate'] = float(np.mean([s[0] for s in stS])); out['gpe_rate'] = float(np.mean([s[0] for s in stG]))
    out['stn_rate_sd'] = float(np.std([s[0] for s in stS]))
    out['stn_cv'] = float(np.nanmean([s[1] for s in stS])); out['gpe_cv'] = float(np.nanmean([s[1] for s in stG]))
    out['stn_burst'] = float(np.mean([s[2] for s in stS])); out['gpe_burst'] = float(np.mean([s[2] for s in stG]))

    edges = np.arange(T0, CFG['T_SIM'] + 1, 1.0)
    for name, sp in (('stn', spS), ('gpe', spG)):
        cnt = np.sum([np.histogram(s, edges)[0] for s in sp], axis=0).astype(float)
        f, P = welch(cnt - cnt.mean(), fs=1000.0, nperseg=1000)
        out[f'{name}_pop_psd'] = P[1:101].tolist()
        out[f'{name}_pop_beta'] = band(f, P, 13, 30)
        out[f'{name}_pop_peak'] = float(f[1:101][np.argmax(P[1:101])])
        out[f'{name}_pop_peak_beta'] = float(f[13:31][np.argmax(P[13:31])])

    if lfp:
        fs = 1000.0 / CFG['REC_DT']
        tI = np.array(ti); mask = tI > T0
        n_t = len(tI)
        I = np.empty((len(imem), int(mask.sum())), dtype=np.float64)
        for k, v in enumerate(imem):
            I[k] = np.asarray(v)[:n_t][mask]
        net = I.sum(0); ab = np.abs(I).sum(0)
        # IClamp current is not part of i_membrane_: it leaves a CONSTANT leftover in the sum, which only
        # shifts the LFP DC level. So the conservation check uses the time-varying part of the sum, and the
        # constant part is compared with the injected bias separately.
        out['cons_ratio'] = float(np.max(np.abs(net - net.mean()) / (ab + 1e-30)))
        out['cons_dc'] = float(net.mean()); out['bias_sum'] = float(sum(ic.amp for ic in ics_stn))
        R = ori
        for ename, E in (('near', CFG['ELEC_NEAR']), ('far', CFG['ELEC_FAR'])):
            E = np.array(E)
            K = np.zeros(len(rows))
            for k, (ci, a, b, r) in enumerate(rows):
                K[k] = line_kernel(R[ci] @ a + pos[ci], R[ci] @ b + pos[ci], r, E, CFG['SIGMA'])
            lf = K @ I
            f, P = psd(lf, fs)
            out[f'{ename}_rms'] = float(lf.std())
            out[f'{ename}_psd'] = P[1:101].tolist()
            for bn, (lo, hi) in BANDS.items():
                out[f'{ename}_{bn}'] = band(f, P, lo, hi)
            if ename == 'near':
                shares = []
                for (a, b) in cell_slices:
                    lc = K[a:b] @ I[a:b]
                    shares.append(band(*psd(lc, fs), 13, 30))
                shares = np.array(shares)
                out['dom_cell_beta_share'] = float(shares.max() / shares.sum())
    del fih
    return out
