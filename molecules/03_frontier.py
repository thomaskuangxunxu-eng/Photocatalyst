"""Step 3: frontier-orbital energies for every catalogue molecule + energy-level diagram.

Usage (from repo root):
    python molecules/03_frontier.py                     # all molecules (~40 min on 10 threads; 4CzIPN ~15 min)
    python molecules/03_frontier.py --only benzophenone DCA
    python molecules/03_frontier.py --plot-only         # redraw the figure from the CSV

Geometry: molecules/xyz_opt/<slug>.xyz if present (step 2), else the MMFF geometry from step 1.
Level: B3LYP/def2-SVP (def2 ECP on I), density fitting, ddCOSMO water, restricted closed shell.
Writes results/molecules/frontier_orbitals.csv and results/molecules/figures/homo_lumo_diagram.png.
Caveats:
  * a Kohn-Sham HOMO-LUMO gap is not an excitation energy (use TDDFT for that);
  * KS orbital energies are not redox potentials. The columns minus_homo_minus_4p44_qualitative /
    minus_lumo_minus_4p44_qualitative are a rough ordering index only (no solvation free energy,
    no relaxation of the ion). For the dye dianions they are off by >1 V from measured potentials.
    For that reason the diagram does not draw the H+/H2 or O2/H2O levels next to the orbital energies.
"""
import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.common import ROOT, build_mol, make_rks, homo_lumo, set_style  # noqa: E402
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--only", nargs="*"); ap.add_argument("--plot-only", action="store_true")
ap.add_argument("--xc", default="b3lyp"); ap.add_argument("--basis", default="def2-svp")
ap.add_argument("--threads", type=int, default=10)
a = ap.parse_args()

cat = pd.read_csv(ROOT / "data/organic_photocatalysts.csv")
out_csv = ROOT / "results/molecules/frontier_orbitals.csv"

if not a.plot_only:
    from pyscf import lib
    lib.num_threads(a.threads)
    todo = cat if not a.only else cat[cat.slug.isin(a.only)]
    todo = todo.sort_values("n_atoms_with_H")
    old = pd.read_csv(out_csv) if out_csv.exists() else pd.DataFrame(columns=["slug"])
    for _, r in todo.iterrows():
        opt = ROOT / f"molecules/xyz_opt/{r.slug}.xyz"
        src = opt if opt.exists() else ROOT / f"molecules/xyz/{r.slug}.xyz"
        t0 = time.time()
        mol = build_mol(src, basis=a.basis)
        mf = make_rks(mol, a.xc); e = mf.kernel()
        if not mf.converged:
            mf = mf.newton(); e = mf.kernel()
        homo, lumo = homo_lumo(mf, mol)
        dip = float(np.linalg.norm(mf.dip_moment(unit="Debye", verbose=0))) if int(r.charge) == 0 else np.nan
        gap = lumo - homo
        row = dict(slug=r.slug, abbrev=r.abbrev, name=r["name"], family=r.family, charge=int(r.charge),
                   geom_source=src.parent.name, nao=mol.nao, n_electrons=mol.nelectron, converged=bool(mf.converged),
                   E_tot_Ha=round(float(e), 6), homo_eV=round(homo, 4), lumo_eV=round(lumo, 4), gap_eV=round(gap, 4),
                   onset_nm=round(1239.841984 / gap, 1), dipole_D=(round(dip, 3) if dip == dip else None),
                   lambda_max_nm=r.lambda_max_nm, lowest_abs_band_nm=r.lowest_abs_band_nm,
                   lowest_abs_band_eV=round(1239.841984 / r.lowest_abs_band_nm, 3),
                   gap_minus_lowest_abs_eV=round(gap - 1239.841984 / r.lowest_abs_band_nm, 3),
                   minus_homo_minus_4p44_qualitative=round(-homo - 4.44, 2),   # ordering index only, NOT V vs SHE
                   minus_lumo_minus_4p44_qualitative=round(-lumo - 4.44, 2),   # ordering index only, NOT V vs SHE
                   wall_s=round(time.time() - t0, 1))
        old = pd.concat([old[old.slug != r.slug], pd.DataFrame([row])], ignore_index=True)
        old.to_csv(out_csv, index=False)
        print(f"{r.slug:14s} nao={mol.nao:5d} conv={mf.converged} HOMO={homo:7.3f} LUMO={lumo:7.3f} "
              f"gap={gap:6.3f} eV  {time.time()-t0:6.1f}s", flush=True)

# ---------------- figure ----------------
import matplotlib.pyplot as plt
set_style()
fo = pd.read_csv(out_csv)
if len(fo) < 2:
    sys.exit(0)
r_ = np.corrcoef(fo.gap_eV, fo.lowest_abs_band_eV)[0, 1]; mae = fo.gap_minus_lowest_abs_eV.abs().mean()
o = fo.sort_values("gap_eV").reset_index(drop=True)
CH = {-2: "#7b3294", 0: "#444444", 1: "#1b7837"}
fig = plt.figure(figsize=(7.6, 4.0))
ax = fig.add_axes([0.085, 0.24, 0.57, 0.67]); ax2 = fig.add_axes([0.745, 0.24, 0.21, 0.5])
for i, rw in o.iterrows():
    c = CH.get(int(rw.charge), "#444444")
    ax.hlines([rw.homo_eV, rw.lumo_eV], i - 0.32, i + 0.32, color=c, lw=2.2)
    ax.vlines(i, rw.homo_eV, rw.lumo_eV, color=c, lw=0.6, ls=(0, (2, 2)))
    ax.text(i, rw.lumo_eV + 0.12, f"{rw.gap_eV:.2f}", ha="center", va="bottom", fontsize=6, color=c)
# No H+/H2 or O2/H2O lines: KS orbital energies are not redox potentials (see docstring).
mmff = o.geom_source == "xyz"
ax.set_xticks(range(len(o)))
ax.set_xticklabels([f"{ab}\u2020" if m else ab for ab, m in zip(o.abbrev, mmff)], rotation=60, ha="right")
ax.set_xlim(-0.6, len(o) - 0.4); ax.set_ylim(-7.7, -0.3)
ax.set_ylabel("KS orbital energy (eV)\nnot a redox potential")
fig.text(0.085, 0.015, f"\u2020 MMFF94s geometry ({int(mmff.sum())} of {len(o)}); "
         f"others B3LYP/def2-SVP optimised (molecules/02_optimize.py)", fontsize=6.5, color="#555", va="bottom")
ax.set_title(f"{a.xc.upper()}/{a.basis.replace('svp', 'SVP').replace('tzvp', 'TZVP')} frontier orbitals (water, ddCOSMO)", loc="left")
for q, lab in [(-2, "dianion"), (0, "neutral"), (1, "cation")]:
    ax.plot([], [], color=CH[q], lw=2.2, label=lab)
ax.legend(loc="upper left", ncol=3)
lim = (1.5, 5.2)
ax2.plot(lim, lim, color="#9a9a9a", lw=0.9, ls="--")
ax2.text(4.05, 4.65, "parity", rotation=45, fontsize=6.5, color="#888", ha="center", va="center")
ax2.scatter(fo.lowest_abs_band_eV, fo.gap_eV, s=22, c=[CH.get(int(q), "#444") for q in fo.charge],
            edgecolor="white", lw=0.5, zorder=3)
for ab, (dx, dy) in {"BP": (-4, 4), "TPP": (5, 4), "MB+": (-4, 7), "Me-PDI": (5, -4)}.items():
    rr = fo[fo.abbrev == ab]
    if len(rr):
        ax2.annotate(ab, (rr.lowest_abs_band_eV.item(), rr.gap_eV.item()), xytext=(dx, dy),
                     textcoords="offset points", fontsize=6.5, ha="left" if dx > 0 else "right")
ax2.set_xlim(lim); ax2.set_ylim(lim); ax2.set_aspect("equal")
ax2.set_xlabel("Lowest abs. band (eV)"); ax2.set_ylabel("KS HOMO\u2013LUMO gap (eV)")
ax2.set_title(f"Gap tracks lowest absorption\nr = {r_:.2f}, MAE {mae:.2f} eV\n(n = {len(fo)})", loc="left")
fig.savefig(ROOT / "results/molecules/figures/homo_lumo_diagram.png", dpi=300)
print(f"figure written; r={r_:.2f}, MAE={mae:.2f} eV, n={len(fo)}")
