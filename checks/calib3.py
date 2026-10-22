import sys, time, json
sys.path.insert(0, '..')
import core4 as c
# half-size network, same in-degree (8 STN<-GPe, 8 GPe<-STN) -> same mean drive per cell
c.CFG.update(N_STN=10, N_GPE=10, P_CONN=0.8, T_SIM=2500.0, T_DISCARD=500.0)
base = dict(ctx_mode='none', bg_rate=100.0, bg_w=0.0025, stn_bias_sd=0.001, gpe_bias_sd=0.002,
            w_gs=0.0004, w_sg=0.004, w_gg=0.0006, d_sg=2.8, d_gs=1.3, d_gg=1.0)
for over in json.loads(sys.argv[1]):
    t = time.time(); kw = dict(base); kw.update(over)
    r = c.simulate(1, lfp=False, **kw)
    print(json.dumps(over), '-> STN %.1f Hz CV %.2f | GPe %.1f Hz CV %.2f | STN pop peak %.0f [%.0fs]' % (
        r['stn_rate'], r['stn_cv'], r['gpe_rate'], r['gpe_cv'], r['stn_pop_peak'], time.time() - t), flush=True)
