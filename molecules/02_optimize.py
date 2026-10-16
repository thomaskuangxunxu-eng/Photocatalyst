"""Step 2: DFT geometry optimization of one molecule with PySCF + geomeTRIC.

Usage (from repo root):
    python molecules/02_optimize.py molecules/xyz/benzene.xyz
    python molecules/02_optimize.py molecules/xyz/eosin_Y.xyz --basis def2-svp --xc b3lyp --disp d3bj --solvent water

Charge is read from the xyz comment line ("charge=-2"), written by 01_represent.py.
Outputs in results/opt/<slug>/: <slug>_opt.xyz, summary.json (energy, HOMO/LUMO, timing).
HOMO/LUMO here are at the optimization level; refine them in 03_frontier.py.
"""
import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
from pyscf import dft, gto, lib
from pyscf.geomopt.geometric_solver import optimize

HARTREE_TO_EV = 27.211386
ROOT = Path(__file__).resolve().parents[1]
HEAVY_ECP = {"I", "Ru", "Ir", "W", "Rh", "Pd", "Pt", "Au", "Ag"}  # def2 ECP elements we may meet


def read_xyz(path: Path) -> tuple[str, int, int]:
    lines = path.read_text().splitlines()
    n = int(lines[0])
    comment = lines[1]
    charge = int(m.group(1)) if (m := re.search(r"charge=(-?\d+)", comment)) else 0
    spin = int(m.group(1)) if (m := re.search(r"spin=(\d+)", comment)) else 0
    atoms = "\n".join(lines[2 : 2 + n])
    return atoms, charge, spin


def build_mol(atoms: str, charge: int, spin: int, basis: str) -> gto.Mole:
    symbols = {line.split()[0] for line in atoms.splitlines()}
    ecp = {el: basis for el in symbols & HEAVY_ECP}  # def2 ECPs share the basis name
    return gto.M(atom=atoms, charge=charge, spin=spin, basis=basis, ecp=ecp or None, verbose=3)


def make_mf(mol: gto.Mole, xc: str, disp: str | None, solvent: str | None):
    mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
    mf.xc = xc
    if disp:
        mf.disp = disp
    if solvent:
        mf = mf.PCM()
        mf.with_solvent.method = "IEF-PCM"
        mf.with_solvent.solvent = solvent
    return mf


def frontier(mf) -> tuple[float, float]:
    e, occ = np.asarray(mf.mo_energy), np.asarray(mf.mo_occ)
    if e.ndim == 2:  # unrestricted: take alpha channel
        e, occ = e[0], occ[0]
    return float(e[occ > 0].max() * HARTREE_TO_EV), float(e[occ == 0].min() * HARTREE_TO_EV)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("xyz", type=Path)
    p.add_argument("--basis", default="def2-svp")
    p.add_argument("--xc", default="b3lyp")
    p.add_argument("--disp", default="d3bj", help="'none' to switch off")
    p.add_argument("--solvent", default=None, help="e.g. water, acetonitrile")
    p.add_argument("--maxsteps", type=int, default=100)
    a = p.parse_args()
    disp = None if a.disp == "none" else a.disp

    slug = a.xyz.stem
    out = ROOT / "results" / "opt" / slug
    out.mkdir(parents=True, exist_ok=True)
    atoms, charge, spin = read_xyz(a.xyz)
    mol = build_mol(atoms, charge, spin, a.basis)
    mol.output = str(out / "pyscf.log")
    mol.build()

    t0 = time.time()
    mf = make_mf(mol, a.xc, disp, a.solvent)
    conv_params = {"convergence_energy": 1e-6, "convergence_grms": 3e-4, "convergence_gmax": 4.5e-4,
                   "convergence_drms": 1.2e-3, "convergence_dmax": 1.8e-3}
    mol_opt = optimize(mf, maxsteps=a.maxsteps, **conv_params)
    mf_final = make_mf(mol_opt, a.xc, disp, a.solvent)
    e_tot = mf_final.kernel()
    wall = time.time() - t0
    homo, lumo = frontier(mf_final)

    coords = mol_opt.atom_coords(unit="Angstrom")
    xyz_lines = [str(mol_opt.natm), f"{slug} {a.xc}/{a.basis} disp={disp} solvent={a.solvent} "
                 f"charge={charge} spin={spin} E={e_tot:.8f} Eh"]
    xyz_lines += [f"{mol_opt.atom_symbol(i):<2} {x:14.8f} {y:14.8f} {z:14.8f}"
                  for i, (x, y, z) in enumerate(coords)]
    (out / f"{slug}_opt.xyz").write_text("\n".join(xyz_lines) + "\n")

    summary = {"name": slug, "xc": a.xc, "basis": a.basis, "disp": disp, "solvent": a.solvent,
               "charge": charge, "spin": spin, "n_atoms": mol_opt.natm, "n_ao": mol_opt.nao,
               "converged": bool(mf_final.converged), "E_tot_Eh": e_tot,
               "HOMO_eV": homo, "LUMO_eV": lumo, "gap_eV": lumo - homo,
               "wall_s": round(wall, 1), "threads": lib.num_threads()}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
