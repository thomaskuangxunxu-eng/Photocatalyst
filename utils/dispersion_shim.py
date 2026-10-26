"""Drop-in replacement for the `pyscf.dispersion` package, backed by the conda-forge
`dftd3-python` (simple-dftd3) and `dftd4-python` libraries.

Why: pyscf-dispersion 1.0.0 (the only release on PyPI) fails on Apple-Silicon Macs for two reasons.
  1. pyscf/dispersion/dftd3.py calls numpy.ctypeslib.load_library('libs-dftd3') without the
     required loader_path argument (TypeError on import).
  2. The wheel ships x86_64 builds of libs-dftd3 / libdftd4, which an arm64 Python cannot load.
pyscf-dispersion is only a thin ctypes wrapper around these same two libraries (Grimme group's
simple-dftd3 and dftd4), so this module re-implements the two classes PySCF uses -
DFTD3Dispersion and DFTD4Dispersion - with the same constructor and get_dispersion(grad) contract,
on top of the official Python bindings.

install() must run BEFORE `import pyscf` (utils/common.py does this).
"""
import sys
import types

import numpy as np


def _numbers_and_coords(mol):
    from pyscf import gto
    numbers = np.asarray([gto.charge(mol.atom_symbol(i)) for i in range(mol.natm)], dtype=np.int32)
    coords = np.asarray(mol.atom_coords(unit="Bohr"), dtype=np.double)   # both libraries expect bohr
    return numbers, coords


def _pack(res, grad):
    out = {"energy": np.array(res["energy"], dtype=np.double)}
    if grad:
        out["gradient"] = np.asarray(res["gradient"], dtype=np.double)
        out["virial"] = np.asarray(res["virial"], dtype=np.double)
    return out


class DFTD3Dispersion:
    """Same signature as pyscf.dispersion.dftd3.DFTD3Dispersion."""
    def __init__(self, mol, xc, version="d3bj", atm=False):
        from dftd3 import interface as d3
        params = {"d3bj": d3.RationalDampingParam, "d3zero": d3.ZeroDampingParam,
                  "d3bjm": d3.ModifiedRationalDampingParam, "d3zerom": d3.ModifiedZeroDampingParam,
                  "d3op": d3.OptimizedPowerDampingParam}
        version = version.lower()
        if version not in params:
            raise ValueError(f"D3 damping '{version}' not supported (use one of {sorted(params)})")
        numbers, coords = _numbers_and_coords(mol)
        self.natm = mol.natm
        self._model = d3.DispersionModel(numbers, coords)
        self._param = params[version](method=xc, atm=bool(atm))

    def get_dispersion(self, grad=False):
        return _pack(self._model.get_dispersion(self._param, grad=grad), grad)


class DFTD4Dispersion:
    """Same signature as pyscf.dispersion.dftd4.DFTD4Dispersion."""
    def __init__(self, mol, xc, atm=False):
        from dftd4 import interface as d4
        numbers, coords = _numbers_and_coords(mol)
        self.natm = mol.natm
        self._model = d4.DispersionModel(numbers, coords, charge=float(mol.charge))
        self._param = d4.DampingParam(method=xc, atm=bool(atm))

    def get_dispersion(self, grad=False):
        return _pack(self._model.get_dispersion(self._param, grad=grad), grad)


def install():
    """Register this module as pyscf.dispersion (+ .dftd3 / .dftd4). Returns True if the
    backing libraries are importable, False otherwise (then nothing is registered)."""
    try:
        import dftd3.interface  # noqa: F401
        import dftd4.interface  # noqa: F401
    except ImportError:
        return False
    pkg = types.ModuleType("pyscf.dispersion")
    m3 = types.ModuleType("pyscf.dispersion.dftd3"); m3.DFTD3Dispersion = DFTD3Dispersion
    m4 = types.ModuleType("pyscf.dispersion.dftd4"); m4.DFTD4Dispersion = DFTD4Dispersion
    pkg.dftd3, pkg.dftd4, pkg.__path__ = m3, m4, []
    pkg.__doc__ = "simple-dftd3 / dftd4 backed replacement (photocatalyst-dft utils/dispersion_shim.py)"
    sys.modules.update({"pyscf.dispersion": pkg, "pyscf.dispersion.dftd3": m3, "pyscf.dispersion.dftd4": m4})
    return True
