"""Experiment definitions for the v4 (published STN cell) network.

States are matched to measured rates from Cruz et al. 2011 (J Neurophysiol 106:2012-2023):
control STN 13.5 Hz / GPe 32.0 Hz; 6-OHDA STN 30.9 Hz / GPe 17.4 Hz.
Only these levers differ between the states (all are assumptions in size, taken from the
direction described by Terman et al. 2002): stronger striatal-like inhibition of GPe
(gpe_bias), weaker GPe-GPe inhibition (w_gg), stronger STN drive (stn_bias).
Delays: STN->GPe 2.8 ms, GPe->STN 1.3 ms, GPe->GPe 1.0 ms (model values of Nevado-Holgado et al. 2014).
"""
import itertools

DELAYS = dict(d_sg=2.8, d_gs=1.3, d_gg=1.0)
COMMON = dict(stn_bias_sd=0.001, gpe_bias_sd=0.002)
NORMAL = dict(w_sg=0.0001, w_gs=0.0004, w_gg=0.0006, gpe_bias=-0.0075, stn_bias=0.004)
PD     = dict(w_sg=0.0001, w_gs=0.0004, w_gg=0.0003, gpe_bias=-0.0185, stn_bias=0.0235)
NOISE  = dict(ctx_mode='none', bg_rate=100.0, bg_w=0.0025)            # non-cortical background only
CORTEX = dict(ctx_mode='poisson', ctx_base=100.0, ctx_w=0.0025, bg_rate=0.0)  # same mean drive, optional beta depth

def _mk(state, extra):
    o = dict(COMMON); o.update(DELAYS); o.update(state); o.update(extra); return o

def exp_B(seeds=tuple(range(1, 11))):
    """State (normal / PD-like) x cortical beta depth (flat, 0.5, 1.0)."""
    J = []
    for s in seeds:
        for name, st in (('N', NORMAL), ('P', PD)):
            for amp in (0.0, 0.5, 1.0):
                J.append(dict(id=f'{name}_amp{amp}_s{s}', seed=s, over=_mk(st, dict(CORTEX, ctx_amp=amp, ctx_freq=20.0))))
    return J

def exp_A(seeds=(1, 2)):
    """Loop only (no cortex), normal state: delay x coupling scale k (w_sg and w_gs both x k)."""
    J = []
    for s in seeds:
        for d, k in itertools.product((3, 5, 10, 15), (1, 3, 6)):
            st = dict(NORMAL); st['w_sg'] *= k; st['w_gs'] *= k
            over = _mk(st, NOISE); over.update(d_sg=float(d), d_gs=float(d))
            J.append(dict(id=f'd{d}_k{k}_s{s}', seed=s, over=over))
    return J

def exp_C(seeds=(1, 2)):
    """What produces the ~50 Hz peak? Change GPe rate, remove one limb of the loop, or both."""
    V = {'base': {}, 'gpe_up': dict(gpe_bias=-0.004), 'gpe_down': dict(gpe_bias=-0.011),
         'no_STN_to_GPe': dict(w_sg=0.0), 'no_GPe_to_STN': dict(w_gs=0.0), 'no_loop': dict(w_sg=0.0, w_gs=0.0)}
    J = []
    for s in seeds:
        for name, ex in V.items():
            st = dict(NORMAL); st.update(ex)
            J.append(dict(id=f'{name}_s{s}', seed=s, over=_mk(st, NOISE)))
    return J

def exp_D(seeds=(1, 2, 3, 4, 5)):
    """Is the ~35 Hz bump of the PD-like state (flat cortex) a loop rhythm or just single-cell firing?
    Same seeds and flat cortex as exp B; remove one or both limbs of the loop."""
    V = {'PD_no_STN_to_GPe': dict(w_sg=0.0), 'PD_no_GPe_to_STN': dict(w_gs=0.0), 'PD_no_loop': dict(w_sg=0.0, w_gs=0.0)}
    J = []
    for s in seeds:
        for name, ex in V.items():
            st = dict(PD); st.update(ex)
            J.append(dict(id=f'{name}_s{s}', seed=s, over=_mk(st, dict(CORTEX, ctx_amp=0.0))))
    return J

def exp_E(seeds=tuple(range(1, 11))):
    """Rate vs network: cross the DRIVE (stn_bias, gpe_bias) with the COUPLING (w_gg) of the two states.
    The two diagonal cells (N drive+N coupling, P drive+P coupling) are experiment B. Here: the off-diagonal cells,
    flat cortex and depth 1.0, same seeds as B."""
    V = {'Pdrive_Ncoup': dict(NORMAL, stn_bias=PD['stn_bias'], gpe_bias=PD['gpe_bias']),
         'Ndrive_Pcoup': dict(PD, stn_bias=NORMAL['stn_bias'], gpe_bias=NORMAL['gpe_bias'])}
    J = []
    for s in seeds:
        for name, st in V.items():
            for amp in (0.0, 1.0):
                J.append(dict(id=f'{name}_amp{amp}_s{s}', seed=s, over=_mk(st, dict(CORTEX, ctx_amp=amp, ctx_freq=20.0))))
    return J

ALL = dict(A=exp_A, B=exp_B, C=exp_C, D=exp_D, E=exp_E)
