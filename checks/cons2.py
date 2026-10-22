import sys, json
sys.path.insert(0,'..')
import core4 as c
c.CFG.update(N_STN=2, N_GPE=2, T_SIM=1200.0, T_DISCARD=600.0)
base = dict(ctx_mode='none', w_gs=0.0, w_sg=0.0, w_gg=0.0)
tests = {'nothing': {}, 'bias_sd': dict(stn_bias_sd=0.001), 'bg': dict(bg_rate=100.0, bg_w=0.0025), 'loop': dict(w_gs=0.0002, w_sg=0.004)}
for k, o in tests.items():
    kw = dict(base); kw.update(o)
    r = c.simulate(1, **kw)
    print(k, 'cons', r['cons_ratio'], 'stn', r['stn_rate'], flush=True)
