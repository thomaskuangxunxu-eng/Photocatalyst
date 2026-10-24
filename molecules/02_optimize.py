"""Step 2: DFT geometry optimisation (PySCF + geomeTRIC) for selected molecules.

Usage (from repo root):
    python molecules/02_optimize.py benzophenone anthraquinone thioxanthone
    python molecules/02_optimize.py eosinY --xc b3lyp --basis def2-svp --solvent water --threads 10
    python molecules/02_optimize.py riboflavin --criteria flat          # flexible molecules, see --help
    python molecules/02_optimize.py eosinY --criteria flat --start partial_last_frame.xyz

Reads molecules/xyz/<slug>.xyz (MMFF geometry from step 1, charge in the comment line).
Writes molecules/xyz_opt/<slug>.xyz and updates results/molecules/geomopt_summary.csv with
cycles, energies and Kabsch-aligned RMSD vs the MMFF geometry.
Default level: B3LYP/def2-SVP, density fitting, ddCOSMO water. D3 is off by default because
pyscf-dispersion is broken in the photocat env (pass --disp d3bj where it works).
"""
import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.common import ROOT, read_xyz, write_xyz, build_mol, make_rks  # noqa: E402  (imports pyscf safely)
import numpy as np
import pandas as pd
from pyscf import lib
from pyscf.geomopt import geometric_solver

BOHR = 0.52917721092


def kabsch_rmsd(P, Q, mask=None):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(P.T @ Q); d = np.sign(np.linalg.det(V @ Wt))
    D = np.sqrt(((P @ (V @ np.diag([1, 1, d]) @ Wt) - Q) ** 2).sum(1))
    if mask is not None: D = D[mask]
    return float(np.sqrt((D ** 2).mean())), float(D.max())


ap = argparse.ArgumentParser()
ap.add_argument("slugs", nargs="+")
ap.add_argument("--xc", default="b3lyp"); ap.add_argument("--basis", default="def2-svp")
ap.add_argument("--solvent", default="water", choices=["water", "none"])
ap.add_argument("--disp", default=None); ap.add_argument("--maxsteps", type=int, default=60)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--criteria", default="gau", choices=["gau", "flat"],
                help="gau: dE 1e-6 Ha, grad rms/max 3e-4/4.5e-4, disp rms/max 1.2e-3/1.8e-3 (geomeTRIC default). "
                     "flat: same gradient thresholds, dE 1e-5 Ha, no displacement thresholds - for molecules with very "
                     "flat torsions (benzoate twist, ribityl chain) where the gau displacement test never settles.")
ap.add_argument("--start", default=None,
                help="xyz to start from (e.g. the last frame of an interrupted run); charge is taken from the MMFF xyz. "
                     "Single slug only. RMSD is always measured against the MMFF geometry.")
a = ap.parse_args()
if a.start and len(a.slugs) != 1:
    sys.exit("--start works with one slug at a time")
CRIT = {"gau": dict(convergence_energy=1e-6, convergence_grms=3e-4, convergence_gmax=4.5e-4,
                    convergence_drms=1.2e-3, convergence_dmax=1.8e-3),
        "flat": dict(convergence_energy=1e-5, convergence_grms=3e-4, convergence_gmax=4.5e-4,
                     convergence_drms=1.0, convergence_dmax=1.0)}[a.criteria]
lib.num_threads(a.threads)
(ROOT / "molecules/xyz_opt").mkdir(parents=True, exist_ok=True)
out_csv = ROOT / "results/molecules/geomopt_summary.csv"
old = pd.read_csv(out_csv) if out_csv.exists() else pd.DataFrame()

for slug in a.slugs:
    t0 = time.time()
    mmff_xyz = ROOT / f"molecules/xyz/{slug}.xyz"
    syms, _, q = read_xyz(mmff_xyz)
    mol0 = build_mol(mmff_xyz, basis=a.basis)                        # MMFF reference geometry
    e0 = make_rks(mol0, a.xc, a.solvent, a.disp).kernel()
    if a.start:
        s_syms, s_xyz, _ = read_xyz(a.start)
        assert s_syms == syms, "start geometry has a different atom order"
        mol = mol0.copy(); mol.set_geom_(np.array(s_xyz), unit="Angstrom")
    else:
        mol = mol0
    mf = make_rks(mol, a.xc, a.solvent, a.disp)
    ncyc = {"n": 0}
    conv, mol_eq = geometric_solver.kernel(mf, maxsteps=a.maxsteps, assert_convergence=False,
                                           callback=lambda env: ncyc.__setitem__("n", ncyc["n"] + 1), **CRIT)
    mf2 = make_rks(mol_eq, a.xc, a.solvent, a.disp); e1 = mf2.kernel()
    c0, c1 = mol0.atom_coords() * BOHR, mol_eq.atom_coords() * BOHR
    heavy = np.array([s != "H" for s in syms])
    rmsd_all, maxdev = kabsch_rmsd(c0, c1); rmsd_heavy, _ = kabsch_rmsd(c0, c1, heavy)
    level = f"{a.xc.upper()}{'-'+a.disp.upper() if a.disp else ''}/{a.basis} {'ddCOSMO(water)' if a.solvent=='water' else 'gas'}"
    write_xyz(ROOT / f"molecules/xyz_opt/{slug}.xyz", syms, c1,
              f"{slug} | charge={q} mult=1 | {level} optimized ({a.criteria} criteria, converged={conv}) | E={e1:.8f} Ha")
    row = dict(slug=slug, level=level, criteria=a.criteria, converged=bool(conv),
               started_from=("MMFF" if not a.start else Path(a.start).name), nao=mol.nao, n_opt_cycles=ncyc["n"],
               E_mmff_geom_Ha=round(e0, 6), E_opt_Ha=round(e1, 6), dE_kcal=round((e1 - e0) * 627.5095, 2),
               wall_s=round(time.time() - t0, 1),
               rmsd_all_A=round(rmsd_all, 4), rmsd_heavy_A=round(rmsd_heavy, 4), max_dev_A=round(maxdev, 4))
    print(row, flush=True)
    old = pd.concat([old[old.slug != slug] if len(old) else old, pd.DataFrame([row])], ignore_index=True)
    old.to_csv(out_csv, index=False)
