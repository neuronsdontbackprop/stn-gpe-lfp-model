import sys, time, json
sys.path.insert(0,'..')
import core4 as c
c.CFG.update(N_STN=2, N_GPE=20, T_SIM=2500.0, T_DISCARD=500.0)
for gb in (0.0, -0.005, -0.01, -0.02, -0.03):
    for wgg in (0.0006, 0.0003):
        r = c.simulate(1, lfp=False, ctx_mode='none', bg_rate=100.0, bg_w=0.0025, w_gs=0.0, w_sg=0.0, w_gg=wgg, gpe_bias=gb, gpe_bias_sd=0.002)
        print('gpe_bias %.3f w_gg %.4f -> GPe %.1f Hz CV %.2f' % (gb, wgg, r['gpe_rate'], r['gpe_cv']), flush=True)
