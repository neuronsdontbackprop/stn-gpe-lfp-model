import sys, os, json, time
ROOT = os.environ.get('V4_ROOT', os.getcwd()); sys.path.insert(0, ROOT)
import core4, jobs4
core4.CFG.update(T_SIM=3000.0, T_DISCARD=500.0)
exp, shard, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
outdir = f'{ROOT}/results/{exp}'; os.makedirs(outdir, exist_ok=True)
J = jobs4.ALL[exp]()
for i, j in enumerate(J):
    if i % n != shard: continue
    fn = f'{outdir}/{j["id"]}.json'
    if os.path.exists(fn): continue
    t = time.time()
    r = core4.simulate(j['seed'], **j['over']); r['id'] = j['id']
    json.dump(r, open(fn + '.tmp', 'w')); os.replace(fn + '.tmp', fn)
    print(j['id'], round(time.time() - t, 1), 's  STN', round(r['stn_rate'], 1), 'GPe', round(r['gpe_rate'], 1), 'cons', r['cons_ratio'], flush=True)
