# updated_work_v4 — what changed, what was found, what is not done

Open `updated_work_v4.ipynb` (already run). All results are in `results/`; set `RUN_SIMS = True` to recompute missing runs (about 160 s each, 2 workers, 111 runs total ≈ 2.5 h).
Tested with NEURON 9.0.2 on Linux. Mod files in `mods/` must be compiled first (`nrnivmodl`; `mknrndll` on Windows). `mods/` holds the 10 channel files of ModelDB 74298 (Gillies & Willshaw 2006) plus two of ours (`gpehh.mod`, `gpecat.mod`); `SThprotocell.hoc` and `sth-data/` are the paper's files (original README: `README_ModelDB_original.txt`). Cite the paper and ModelDB 74298.

## Done in v4
| Professor's point | What was done |
|---|---|
| Use a published ModelDB cell and name it | STN cells = Gillies & Willshaw 2006, ModelDB 74298, rat STN projection neuron. Compiles, fires on its own (8.9 Hz at 35 °C), rebound burst reproduced (`checks/`). |
| LFP from the network's own STN cells | LFP computed from membrane currents of exactly the 20 simulated published cells (3-D layout built from the cell's own sections). |
| Intrinsic activity, normal vs Parkinson | Pacemaking cell (CV 0.03), network rates and CV reported per state. |
| Realistic inputs / different rates per state | States set from Cruz et al. 2011 (STN 13.5 → 30.9 Hz, GPe 32.0 → 17.4 Hz); model gives 14–15 → 31–32 Hz (STN) and 28 → 15 Hz (GPe). Delays 2.8/1.3 ms from Nevado-Holgado et al. 2014. Both papers' numbers checked by hand. |
| Temperature | 35 °C everywhere (published cell native; GPe rates referenced to 35 °C). |
| Vary coupling and delay | Experiment A (delay 3–15 ms × coupling ×1/×3/×6, cortex off). |
| Cortex must exist; beta must not be tuned in | Experiment B: cortical beta depth 0 / 0.5 / 1.0 at 20 Hz, crossed with state, 10 seeds. Loop alone does not make beta at measured delays. |
| Test the 49 Hz harmonic | Experiments C and D: peak is not at a fixed 49 Hz, survives cutting the loop, tracks the STN firing rate. Source not fully identified. |
| ISI CV / bursting | Reported (STN CV 0.38 normal, 0.23 PD-like; burst index 0). |
| Statistics | Log power, paired two-sided Wilcoxon + sign tests, Holm per electrode, 10 seeds. |
| Conservation check | Redone: tests the time-varying part (bias-electrode current is not in `i_membrane_`); 1e-14 in all B runs. v3 Exp. A value of 0.6–0.94 explained (constant offset only). |

## Main findings (near electrode, Holm-corrected)
- PD-like state, flat cortex: no extra LFP beta (×1.10, p = 0.56).
- Cortical beta depth raises LFP beta in the PD-like state (×1.48, ×1.86; Holm p = 0.016; 10/10 seeds), weakly in normal (×1.09, ×1.25; ns).
- State difference at depth 1.0: ×2.22 (Holm p = 0.049); at depth 0.5: ×1.49 (Holm p = 0.078, ns). Sign tests weaker (8/10, p = 0.11). Interaction ×2.02 (Holm p = 0.035). Far electrode: nothing significant.
- Loop alone: peak frequency follows delay (52/39/26/14 Hz for 3/5/10/15 ms), beta only at ~10 ms each way; stronger coupling silences the STN.
- 30–40 Hz bump in the PD-like state is single-cell firing, not a loop rhythm (survives cutting the loop, moves with firing rate).

## Not done / limits
GPe is not a published model (reduced HH point cell); GPe rates 12–15 % below measured; size of PD-like changes is an assumption and the coupling weights are the same in both states (the states differ in drive and GPe–GPe weight), so rate and amplification of cortical beta are confounded; no increase of bursting in the PD-like state; no STN–STN collaterals or striatal input; cortical beta is a sinusoidal modulation (no bursts); invented dendritic 3-D directions; 20 cells; A and C have 2 seeds; phase-randomised cortex (`ctx_phase_jitter`) not run; GABA-only / AMPA-only receptor tests not run (only limb-removal tests); no comparison with in-vivo LFP values; no requirements.txt; novelty statement and citations (Hagen 2016, LFPy, LFPsim/Parasuram 2016, Gillies & Sterratt) still to be written in the report.

## Added: experiment E (rate vs network coupling)
40 runs (`jobs4.exp_E`, results in `results/E`, analysis in `analysis_E.py`). Total is now 151 runs. The 2x cortical-beta gain of the Parkinson-like state comes from the firing-rate drive (stn_bias, gpe_bias; ratio 2.2x, Holm 0.010 / 0.023), not from the GPe-GPe coupling w_gg (0.92x, Holm 0.55).
