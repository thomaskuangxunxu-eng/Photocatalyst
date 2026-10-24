"""Step 4: orbital cube files, fragment localisation, and isosurface figure.

Usage (from repo root):
    python molecules/04_orbitals.py                  # Mes-Acr+ and eosin Y (default figure), ~4 min
    python molecules/04_orbitals.py --slugs DCA --no-figure
    python molecules/04_orbitals.py --plot-only      # redraw from existing cubes + CSV

For each molecule: B3LYP/def2-SVP + ddCOSMO water on the step-1 (or step-2) geometry, then
  * cube files HOMO-2 .. LUMO+1 -> results/molecules/cubes/<slug>_<orbital>.cube  (git-ignored, ~3 MB each)
  * Mulliken orbital population on a "core" fragment defined by a SMARTS pattern in FRAGMENTS
    -> results/molecules/orbital_localization.csv
The figure (results/molecules/figures/frontier_orbitals_isosurfaces.png) shows Mes-Acr+ HOMO-1 / LUMO
and eosin Y HOMO / LUMO. Rendering uses marching cubes in matplotlib (py3Dmol has no headless PNG export).
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils.common import ROOT, HARTREE_EV, build_mol, make_rks, set_style  # noqa: E402
import numpy as np
import pandas as pd
from rdkit import Chem

# slug -> (SMARTS of the core fragment, core name, remainder name)
FRAGMENTS = {
    "MesAcr": ("c1ccc2[n+](C)c3ccccc3cc2c1", "acridinium", "mesityl"),
    "eosinY": ("c1cc(=O)cc2oc3cc([O-])ccc3cc12", "xanthene core", "benzoate ring"),
    "roseBengal": ("c1cc(=O)cc2oc3cc([O-])ccc3cc12", "xanthene core", "benzoate ring"),
    "fluorescein": ("c1cc(=O)cc2oc3cc([O-])ccc3cc12", "xanthene core", "benzoate ring"),
    "4CzIPN": ("N#Cc1cccc(C#N)c1", "isophthalonitrile", "carbazoles"),
}

ap = argparse.ArgumentParser()
ap.add_argument("--slugs", nargs="*", default=["MesAcr", "eosinY"])
ap.add_argument("--no-figure", action="store_true"); ap.add_argument("--threads", type=int, default=10)
ap.add_argument("--grid", type=int, default=70)
ap.add_argument("--plot-only", action="store_true", help="redraw the figure from existing cubes/CSV")
a = ap.parse_args()

cat = pd.read_csv(ROOT / "data/organic_photocatalysts.csv").set_index("slug")
cube_dir = ROOT / "results/molecules/cubes"; cube_dir.mkdir(parents=True, exist_ok=True)
loc_csv = ROOT / "results/molecules/orbital_localization.csv"
from pyscf import lib
from pyscf.tools import cubegen
lib.num_threads(a.threads)

rows = []
for slug in ([] if a.plot_only else a.slugs):
    opt = ROOT / f"molecules/xyz_opt/{slug}.xyz"
    src = opt if opt.exists() else ROOT / f"molecules/xyz/{slug}.xyz"
    mol = build_mol(src); mf = make_rks(mol); mf.kernel()
    no = mol.nelec[0]
    core, cname, rname = set(), "core", "rest"
    if slug in FRAGMENTS:
        smarts, cname, rname = FRAGMENTS[slug]
        rm = Chem.AddHs(Chem.MolFromSmiles(cat.loc[slug, "smiles"]))
        core = set(rm.GetSubstructMatch(Chem.MolFromSmarts(smarts)))
        if not core:
            sys.exit(f"fragment SMARTS for {slug} does not match its SMILES - fix FRAGMENTS")
        for at in rm.GetAtoms():   # terminal atoms (H, halogens, O) follow the atom they hang on
            if at.GetIdx() not in core and at.GetDegree() == 1 and at.GetNeighbors()[0].GetIdx() in core:
                core.add(at.GetIdx())
    S = mf.get_ovlp(); sl = mol.aoslice_by_atom()
    for k in range(no - 3, no + 2):
        lab = "HOMO" if k == no - 1 else "LUMO" if k == no else (f"HOMO-{no-1-k}" if k < no - 1 else f"LUMO+{k-no}")
        c = mf.mo_coeff[:, k]; q = c * (S @ c)
        frac = sum(q[a0:a1].sum() for i, (_, _, a0, a1) in enumerate(sl) if i in core) if core else np.nan
        cubegen.orbital(mol, str(cube_dir / f"{slug}_{lab}.cube"), c, nx=a.grid, ny=a.grid, nz=a.grid)
        rows.append(dict(slug=slug, orbital=lab, E_eV=round(mf.mo_energy[k] * HARTREE_EV, 3),
                         core_fragment=cname, pct_core=round(100 * frac, 1), other_fragment=rname,
                         pct_other=round(100 - 100 * frac, 1)))
    print(pd.DataFrame([r for r in rows if r["slug"] == slug]).to_string(index=False), flush=True)

old = pd.read_csv(loc_csv) if loc_csv.exists() else pd.DataFrame(columns=["slug"])
new = pd.DataFrame(rows, columns=["slug"]) if not rows else pd.DataFrame(rows)
pd.concat([old[~old.slug.isin(new.slug)], new], ignore_index=True).to_csv(loc_csv, index=False)

if a.no_figure:
    sys.exit(0)

# ---------------- isosurface figure ----------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes
set_style()
B2A = 0.52917721092
ELC = {1: "#d9d9d9", 6: "#4d4d4d", 7: "#3050f8", 8: "#ff2d2d", 35: "#a62929", 53: "#8f00bf"}


def read_cube(p):
    L = Path(p).read_text().split("\n"); nat = int(L[2].split()[0]); org = np.array(L[2].split()[1:4], float)
    axes = [(int(L[3 + i].split()[0]), np.array(L[3 + i].split()[1:4], float)) for i in range(3)]
    atoms = [(int(L[6 + i].split()[0]), np.array(L[6 + i].split()[2:5], float)) for i in range(nat)]
    vals = np.array(" ".join(L[6 + nat:]).split(), float).reshape([n for n, _ in axes])
    return org, [v for _, v in axes], atoms, vals


def render(ax, cube, smiles, title, iso=0.035, elev=50, azim=-75):
    org, axes, atoms, vals = read_cube(cube)
    step = np.array([np.linalg.norm(v) for v in axes])
    P = np.array([p for _, p in atoms]); Zs = np.array([z for z, _ in atoms])
    ctr = P.mean(0); _, _, Vt = np.linalg.svd(P - ctr); tf = lambda X: (X - ctr) @ Vt.T * B2A
    for sgn, col in ((1, "#1b6ca8"), (-1, "#d9822b")):
        verts, faces, *_ = marching_cubes(sgn * vals, iso, spacing=tuple(step))
        ax.add_collection3d(Poly3DCollection(tf(verts + org)[faces], facecolor=col, edgecolor="none", alpha=0.5))
    Q = tf(P); m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); ax.plot(*Q[[i, j]].T, color="#444", lw=1.1)
    hv = Zs > 1; ax.scatter(*Q[hv].T, s=7, c=[ELC.get(z, "#999") for z in Zs[hv]], depthshade=False)
    lo, hi = Q.min(0) - 1.2, Q.max(0) + 1.2; c = (lo + hi) / 2; h = (hi - lo).max() / 2
    ax.set_xlim(c[0] - h, c[0] + h); ax.set_ylim(c[1] - h, c[1] + h); ax.set_zlim(c[2] - h, c[2] + h)
    ax.set_proj_type("ortho"); ax.set_box_aspect((1, 1, 1), zoom=1.45)
    ax.view_init(elev=elev, azim=azim); ax.set_axis_off()
    ax.text2D(0.5, 0.97, title, transform=ax.transAxes, ha="center", va="top", fontsize=8)


loc = pd.read_csv(loc_csv).set_index(["slug", "orbital"])
need = [("MesAcr", "HOMO-1"), ("MesAcr", "LUMO"), ("eosinY", "HOMO"), ("eosinY", "LUMO")]
if not all((s, o) in loc.index for s, o in need):
    print("figure skipped: run with --slugs MesAcr eosinY first"); sys.exit(0)
L = lambda s, o, col: loc.loc[(s, o), col]
gap_ma = 1000 * (L("MesAcr", "HOMO", "E_eV") - L("MesAcr", "HOMO-1", "E_eV"))
panels = [("MesAcr", "HOMO-1", f"Mes-Acr$^+$ HOMO\u22121 ({L('MesAcr','HOMO-1','E_eV'):.2f} eV)\n{L('MesAcr','HOMO-1','pct_other'):.0f}% on mesityl (donor)", 25, -50),
          ("MesAcr", "LUMO", f"Mes-Acr$^+$ LUMO ({L('MesAcr','LUMO','E_eV'):.2f} eV)\n{L('MesAcr','LUMO','pct_core'):.0f}% on acridinium (acceptor)", 25, -50),
          ("eosinY", "HOMO", f"Eosin Y$^{{2-}}$ HOMO ({L('eosinY','HOMO','E_eV'):.2f} eV)\n{L('eosinY','HOMO','pct_core'):.0f}% on xanthene core", 55, -80),
          ("eosinY", "LUMO", f"Eosin Y$^{{2-}}$ LUMO ({L('eosinY','LUMO','E_eV'):.2f} eV)\n{L('eosinY','LUMO','pct_core'):.0f}% on xanthene core", 55, -80)]
fig = plt.figure(figsize=(7.2, 6.8))
for k, (s, o, title, el, az) in enumerate(panels):
    axk = fig.add_subplot(2, 2, k + 1, projection="3d")
    render(axk, cube_dir / f"{s}_{o}.cube", cat.loc[s, "smiles"], title, elev=el, azim=az)
fig.suptitle("Mes-Acr$^+$ frontier orbitals sit on orthogonal rings; eosin Y orbitals share the xanthene", y=0.99, fontsize=9)
fig.text(0.5, 0.012, "Isosurfaces \u00b10.035 a.u. (blue +, orange \u2212); B3LYP/def2-SVP, ddCOSMO water.\n"
         f"The Mes-Acr$^+$ HOMO is acridinium-centred and lies only {gap_ma:.0f} meV above this mesityl orbital.",
         ha="center", fontsize=7)
fig.subplots_adjust(left=0.0, right=1.0, top=0.95, bottom=0.04, wspace=0.0, hspace=0.02)
fig.savefig(ROOT / "results/molecules/figures/frontier_orbitals_isosurfaces.png", dpi=250)
print("figure written")
