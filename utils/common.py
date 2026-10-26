"""Shared helpers for the photocatalyst-dft pipeline.

Import this module BEFORE pyscf in every script:  from utils.common import ...
"""
import importlib
import re
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HARTREE_EV = 27.211386245988
WATER_EPS = 78.3553          # ddCOSMO dielectric constant for water
H2_LEVEL_EV = -4.44           # H+/H2 on the absolute (vacuum) scale
O2_LEVEL_EV = -5.67           # O2/H2O on the absolute scale


def _purge_pyscf():
    for k in [k for k in sys.modules if k == "pyscf" or k.startswith("pyscf.")]:
        del sys.modules[k]


def _import_pyscf():
    """Import pyscf with a working D3/D4 dispersion backend where possible.

    Order of preference:
      1. utils/dispersion_shim.py - pyscf.dispersion re-implemented on the conda-forge
         `dftd3-python` / `dftd4-python` bindings (works on Apple Silicon).
      2. the PyPI `pyscf-dispersion` package, if it imports cleanly (e.g. x86-64 Linux).
      3. no dispersion: an empty pyscf.dispersion placeholder so `import pyscf` still works.
    pyscf-dispersion 1.0.0 is broken on arm64 Macs (load_library() call without loader_path,
    and x86_64-only binaries), which is why the shim comes first.
    Returns (pyscf module, backend name or None).
    """
    if "pyscf" in sys.modules:   # imported before utils.common: cannot swap the backend any more
        try:
            from pyscf.dispersion import dftd3  # noqa: F401
            return sys.modules["pyscf"], "already-imported pyscf.dispersion"
        except Exception:
            return sys.modules["pyscf"], None
    from utils import dispersion_shim
    if dispersion_shim.install():
        return importlib.import_module("pyscf"), "simple-dftd3 / dftd4 (utils/dispersion_shim.py)"
    try:
        mod = importlib.import_module("pyscf")
        from pyscf.dispersion import dftd3  # noqa: F401
        return mod, "pyscf-dispersion"
    except ImportError:
        return importlib.import_module("pyscf"), None      # no dispersion package at all
    except (TypeError, OSError):                            # broken pyscf-dispersion build
        _purge_pyscf()
        sys.modules["pyscf.dispersion"] = types.ModuleType("pyscf.dispersion")
        return importlib.import_module("pyscf"), None


pyscf, DISPERSION_BACKEND = _import_pyscf()
DISPERSION_AVAILABLE = DISPERSION_BACKEND is not None


def read_xyz(path):
    """Return (symbols, coords[Angstrom] as list of tuples, charge) from a .xyz file.
    Charge is parsed from 'charge=<int>' in the comment line (default 0)."""
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    m = re.search(r"charge=(-?\d+)", lines[1])
    charge = int(m.group(1)) if m else 0
    syms, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        syms.append(p[0]); xyz.append(tuple(float(v) for v in p[1:4]))
    return syms, xyz, charge


def write_xyz(path, syms, coords, comment):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{comment}\n")
        for s, (x, y, z) in zip(syms, coords):
            f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def build_mol(xyz_path, basis="def2-svp", max_memory=6000):
    """PySCF Mole from an xyz written by this pipeline; def2 ECPs for Z > 36 (I)."""
    from pyscf import gto
    from pyscf.data.elements import charge as Z
    syms, xyz, q = read_xyz(xyz_path)
    use_ecp = any(Z(s) > 36 for s in syms)
    return gto.M(atom=[(s, c) for s, c in zip(syms, xyz)], basis=basis,
                 ecp=(basis if use_ecp else None), charge=q, spin=0,
                 verbose=0, max_memory=max_memory)


def make_rks(mol, xc="b3lyp", solvent="water", disp=None, conv_tol=1e-8):
    """Density-fitted restricted KS, optionally with ddCOSMO water and D3."""
    from pyscf import dft
    mf = dft.RKS(mol, xc=xc).density_fit()
    if disp:
        if not DISPERSION_AVAILABLE:
            raise RuntimeError("no dispersion backend: `conda install -c conda-forge dftd3-python dftd4-python` "
                               "(see utils/dispersion_shim.py), or run without --disp")
        mf.disp = disp
    if solvent == "water":
        mf = mf.ddCOSMO(); mf.with_solvent.eps = WATER_EPS
    mf.conv_tol = conv_tol; mf.max_cycle = 100
    return mf


def homo_lumo(mf, mol):
    no = mol.nelec[0]
    return mf.mo_energy[no - 1] * HARTREE_EV, mf.mo_energy[no] * HARTREE_EV


def set_style():
    """Matplotlib defaults used for the figures in results/."""
    import matplotlib as mpl
    mpl.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "xtick.direction": "out", "ytick.direction": "out",
                         "legend.frameon": False, "savefig.dpi": 300, "figure.dpi": 100})
