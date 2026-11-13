import sys, time, json, warnings, numpy as np
warnings.filterwarnings("ignore")
import pyDeltaRCM
nsteps = int(sys.argv[1]); seeds = [int(s) for s in sys.argv[2].split(",")]
out = {}
for seed in seeds:
    m = pyDeltaRCM.DeltaModel(out_dir=f"{sys.argv[3]}/out_{seed}", seed=seed)
    m.update()  # warm up (numba compile)
    t0 = time.perf_counter()
    for _ in range(nsteps - 1): m.update()
    dt = time.perf_counter() - t0
    eta = m.eta.astype(float); ei = m.eta_init.astype(float)
    land = (eta > -0.5*m.h0) & (m.cell_type == 0)
    wet = m.depth > m.dry_depth
    rows = np.where(land.any(axis=1))[0] if land.any() else np.array([0])
    out[seed] = dict(sec_per_step=dt/(nsteps-1),
        shore=int(rows.max()), land_cells=int(land.sum()),
        vol=float((eta-ei)[m.cell_type!=-2].sum()*m.dx**2),
        mean_uw=float(m.uw[wet].mean()), mean_depth=float(m.depth[wet].mean()),
        mean_stage_land=float(m.stage[land].mean()) if land.any() else 0.0,
        profile=(np.maximum(eta-ei,0).sum(axis=1)*m.dx**2).tolist(),
        sandfrac_mean=float(m.sand_frac[(eta-ei)>0.1].mean()))
json.dump(out, open(f"{sys.argv[3]}/py_{nsteps}_{sys.argv[2]}.json","w"))
print({k:(round(v['sec_per_step'],3), v['shore'], v['land_cells'], round(v['vol'])) for k,v in out.items()})
