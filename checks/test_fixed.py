import sys, time
from neuron import h
h.load_file("stdgui.hoc")
h.xopen("SThprotocell.hoc")
h('objref SThtype, SThcells[1]')
h.SThtype = h.SThproto()
h.SThcells[0] = h.SThcell(0, h.SThtype)
h.xopen("tools.hoc")
c = h.SThcells[0]
c.soma.gna_Na = 1.483419823e-02
c.soma.gna_NaL = 1.108670852e-05
for nm in ["gk_KDR","gk_Kv31","gk_Ih","gk_sKCa","gcaT_CaT","gcaN_HVA","gcaL_HVA"]:
    h.cset(0, nm, "")
h('forall ion_style("na_ion",1,2,1,0,1)')
h('forall ion_style("k_ion",1,2,1,0,1)')
h('forall ion_style("ca_ion",3,2,1,1,1)')
h.set_aCSF(4)
temp = float(sys.argv[1]) if len(sys.argv)>1 else 37
h.celsius = temp
cv = h.CVode(); cv.active(0)
cv.atolscale("cai", 1.0e-3)
tv = h.Vector(); vv = h.Vector()
tv.record(h._ref_t); vv.record(c.soma(0.5)._ref_v)
apc = h.APCount(c.soma(0.5)); apc.thresh = -20
st = h.Vector(); apc.record(st)
h.tstop = 1000; h.dt = 0.025
t0 = time.time(); h.v_init=-65; h.stdinit(); h.run()
print("celsius", temp, "wall s", round(time.time()-t0,1), "nsec", sum(1 for _ in h.allsec()), "nseg", sum(s.nseg for s in h.allsec()))
sp = [x for x in st if x > 200]
print("spikes after 200ms:", len(sp), "rate Hz", len(sp)/1.9)
import numpy as np
if len(sp) > 3:
    isi = np.diff(sp); print("ISI mean", isi.mean().round(1), "CV", (isi.std()/isi.mean()).round(3))
v = np.array(vv); t = np.array(tv)
print("v min/max after 200ms", v[t>200].min().round(1), v[t>200].max().round(1))
