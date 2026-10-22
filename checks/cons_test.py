import sys
sys.path.insert(0,'..')
import numpy as np
import core4 as c
from neuron import h
c.proto(); h.cvode.use_fast_imem(1)
mode = sys.argv[1]; dt = float(sys.argv[2])
cell, secs = c.build_stn(0)
if mode == 'nostyle':
    for s in secs:
        h.ion_style('na_ion',0,1,0,0,0,sec=s); h.ion_style('k_ion',0,1,0,0,0,sec=s); h.ion_style('ca_ion',0,1,0,0,0,sec=s)
h.nai0_na_ion = 15; h.nao0_na_ion = 128.5; h.ki0_k_ion = 140; h.ko0_k_ion = 2.5; h.cai0_ca_ion=1e-4; h.cao0_ca_ion=2.0
vecs=[h.Vector().record(seg._ref_i_membrane_, 0.25) for s in secs for seg in s]
h.dt=dt; h.celsius=35; h.finitialize(-65); h.continuerun(600)
I=np.array([np.array(v) for v in vecs]); I=I[:,400:]
net=I.sum(0); ab=np.abs(I).sum(0)
print(mode, dt, 'max|net| nA', np.abs(net).max(), 'max sum|I| nA', ab.max(), 'ratio max', (np.abs(net)/ab).max(), 'ratio of maxes', np.abs(net).max()/ab.max())
