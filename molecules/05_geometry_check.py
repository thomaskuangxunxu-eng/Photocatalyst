"""Step 5: how much do the frontier orbitals depend on MMFF vs DFT geometry?

For every molecule that has a DFT geometry in molecules/xyz_opt/, run the step-3 single point
(B3LYP/def2-SVP, ddCOSMO water) on the MMFF geometry from molecules/xyz/ and compare it with the
DFT-geometry values already in results/molecules/frontier_orbitals.csv (run step 3 after step 2).
Writes results/molecules/geometry_sensitivity.csv with the HOMO / LUMO / gap shifts (DFT minus MMFF)
next to the geometry change from results/molecules/geomopt_summary.csv.

Usage (from repo root):  python molecules/05_geometry_check.py [--threads 10]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.common import ROOT, build_mol, make_rks, homo_lumo  # noqa: E402  (imports pyscf safely)
import pandas as pd
from pyscf import lib

ap = argparse.ArgumentParser(); ap.add_argument("--threads", type=int, default=10)
a = ap.parse_args(); lib.num_threads(a.threads)

fo = pd.read_csv(ROOT / "results/molecules/frontier_orbitals.csv").set_index("slug")
gs = pd.read_csv(ROOT / "results/molecules/geomopt_summary.csv")
gs = gs[~gs.level.str.contains(r"-D[34]", regex=True)].set_index("slug")   # xyz_opt/ holds the no-dispersion level
rows = []
for xyz in sorted((ROOT / "molecules/xyz_opt").glob("*.xyz")):
    slug = xyz.stem
    if slug not in fo.index or fo.loc[slug, "geom_source"] != "xyz_opt":
        print(f"skip {slug}: frontier_orbitals.csv not computed on the DFT geometry (rerun step 3)"); continue
    mol = build_mol(ROOT / f"molecules/xyz/{slug}.xyz"); mf = make_rks(mol); mf.kernel()
    h, l = homo_lumo(mf, mol); d = fo.loc[slug]
    rows.append(dict(slug=slug, abbrev=d.abbrev, charge=int(d.charge),
                     rmsd_heavy_A=gs.loc[slug, "rmsd_heavy_A"], max_dev_A=gs.loc[slug, "max_dev_A"],
                     dE_opt_kcal=gs.loc[slug, "dE_kcal"],
                     homo_mmff_eV=round(h, 4), lumo_mmff_eV=round(l, 4), gap_mmff_eV=round(l - h, 4),
                     homo_dft_eV=d.homo_eV, lumo_dft_eV=d.lumo_eV, gap_dft_eV=d.gap_eV,
                     dHOMO_eV=round(d.homo_eV - h, 3), dLUMO_eV=round(d.lumo_eV - l, 3),
                     dGap_eV=round(d.gap_eV - (l - h), 3), scf_converged=bool(mf.converged)))
    print(rows[-1], flush=True)
out = pd.DataFrame(rows)
out.to_csv(ROOT / "results/molecules/geometry_sensitivity.csv", index=False)
print(out[["abbrev", "rmsd_heavy_A", "dE_opt_kcal", "dHOMO_eV", "dLUMO_eV", "dGap_eV"]].to_string(index=False))
