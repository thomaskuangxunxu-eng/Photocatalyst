# photocatalyst-dft

Computational survey of reported photocatalysts: identifiers, 3D structures and DFT frontier orbitals for
organic dyes (RDKit + PySCF). Solid photocatalysts (Materials Project + periodic DFT) are done and will be
added to this repository in a later update.

The results and what they mean are in **[REPORT.md](REPORT.md)**. This README covers how to rerun them.

## Setup

```bash
conda env create -f environment.yml
conda activate photocat
```

All scripts run from the repository root and import `utils/common.py` before PySCF.

**Dispersion (D3/D4).** Do not install `pyscf-dispersion`. Version 1.0.0 calls NumPy's
`load_library()` without its required path argument and ships x86_64-only libraries, so it breaks
`import pyscf` on Apple Silicon. Instead, `utils/dispersion_shim.py` provides `pyscf.dispersion` on top of
the conda-forge `dftd3-python` / `dftd4-python` bindings (the same Grimme-group libraries, native arm64),
and `utils/common.py` loads it automatically. Use `--disp d3bj` (or `d3zero`, `d4`) with
`molecules/02_optimize.py`; D3 geometries go to `molecules/xyz_opt_d3bj/`.
- **Validation:** D3 is additive to the SCF energy to 10⁻¹⁵ Ha, and analytic and finite-difference gradients agree to 4×10⁻¹³ Ha/bohr.
- **No effect on orbitals at a fixed geometry.** D3 is a post-SCF energy term, so the HOMO/LUMO energies there are unchanged.
- **The stored results are without dispersion.** For benzophenone and thioxanthone, re-optimising with
  D3(BJ) moved the geometry by ≤ 0.022 Å heavy-atom RMSD and the gap by ≤ 0.006 eV.
  The flexible dyes (eosin Y, riboflavin, rhodamine B) have not been re-optimised with D3 yet.

## Layout

| Path | Contents |
| --- | --- |
| `data/organic_photocatalysts.csv` | 16 organic photocatalysts: SMILES, InChIKey, PubChem parent CID (+ how it relates to the stored ion), family, role, literature λmax, lowest absorption band, redox notes |
| `molecules/xyz/` | MMFF94s geometries (step 1), charge in the comment line |
| `molecules/xyz_opt/` | B3LYP/def2-SVP + ddCOSMO(water) optimised geometries (step 2) |
| `results/molecules/` | tables (CSV), `figures/`, 3D viewer (`organic_photocatalysts_3d.html`) |
| `utils/common.py` | shared helpers (safe PySCF import, xyz I/O, KS setup, plot style) |

## Pipeline — molecules

| Step | Command | Output | Cost (12-core Mac) |
| --- | --- | --- | --- |
| 1. Represent | `python molecules/01_represent.py` | checks formula/charge/InChIKey; `molecules/xyz/*.xyz`; grid figure; `results/molecules/organic_photocatalysts_3d.html` | seconds |
| 2. Optimise | `python molecules/02_optimize.py benzophenone anthraquinone thioxanthone`<br>`python molecules/02_optimize.py riboflavin --criteria flat` | `molecules/xyz_opt/`, `results/molecules/geomopt_summary.csv` | rigid: 4–7 min; eosin Y ≈ 1.4 h; riboflavin ≈ 1.5 h; rhodamine B ≈ 19 h |
| 3. Frontier orbitals | `python molecules/03_frontier.py [--only slug ...] [--plot-only]` | `results/molecules/frontier_orbitals.csv`, `figures/homo_lumo_diagram.png` | ~40 min all (4CzIPN ~15 min) |
| 4. Orbitals | `python molecules/04_orbitals.py [--slugs ...] [--plot-only]` | cubes (git-ignored), `orbital_localization.csv`, `figures/frontier_orbitals_isosurfaces.png` | ~4 min |
| 5. Geometry check | `python molecules/05_geometry_check.py` | `results/molecules/geometry_sensitivity.csv`: HOMO/LUMO/gap on MMFF vs DFT geometry for every optimised molecule | ~8 min |

Level of theory: B3LYP/def2-SVP (def2 ECP on I), density fitting, ddCOSMO water (ε = 78.36). Steps 3–4 use
the step-2 DFT geometry when one exists, otherwise the MMFF94s geometry. Six molecules have DFT geometries:
TX, BP, AQ, EY2-, RF, RhB+. The other 10 (DCA, BZ, PTH, RB2-, FL2-, MB+, Mes-Acr+, Me-PDI, TPP, 4CzIPN) are still on MMFF geometries and are marked † in the
diagram. Step 5 measures what that costs. For the rigid molecules the gap changes by ≤ 0.03 eV between MMFF and DFT
geometry. For the flexible or charged ones it changes by up to +0.15 eV (EY2-), with heavy-atom RMSD 0.17–0.40 Å.
RB²⁻ and FL²⁻ share eosin Y's scaffold and probably carry a similar error. Use `--criteria flat` for molecules with very
flat torsions, where geomeTRIC's default displacement test never settles (see `python molecules/02_optimize.py --help`).

`figures/structures_3d_montage.png` and `figures/geomopt_overlay_benzophenone.png` were drawn interactively
and have no script.

## Key results (details in REPORT.md)

- **Organic:** all 16 B3LYP SCFs converged. The KS gap tracks the lowest experimental absorption band
  (r = 0.89, MAE 0.46 eV, n = 16). Benzophenone is the outlier (+1.37 eV), because its first band is n→π*.
  The KS gap is not an excitation energy; TDDFT is the next step.
- **Orbital energies are not redox potentials.** The diagram therefore has no H⁺/H₂ or O₂/H₂O lines. The columns
  `minus_homo_minus_4p44_qualitative` / `minus_lumo_minus_4p44_qualitative` in `frontier_orbitals.csv` are an ordering index
  only; for the dianions they are off by more than 1 V from measured potentials.
- **Geometry:** for flexible or charged dyes, MMFF geometries shift the gap by up to 0.15 eV (`geometry_sensitivity.csv`).
- **Mes-Acr⁺:** the LUMO is 98 % acridinium. The mesityl donor orbitals (HOMO−1/−2) lie within 13 meV of an
  acridinium-centred HOMO, so the KS ordering of those three is not meaningful at this level.

## Sanity checks

| Check | Result | Status |
| --- | --- | --- |
| 16 SMILES parse; formula, charge, InChIKey match catalogue | 16/16 (`01_represent.py` exits non-zero otherwise) | pass |
| `01_represent.py` reproducible | XYZ byte-identical on rerun (fixed ETKDG seed) | pass |
| `03_frontier.py` rerun (BP, TX, DCA) vs stored values | max difference 0.000 eV | pass |
| `02_optimize.py` rerun (thioxanthone) | same energy (−973.119144 Ha), 10 cycles, heavy-atom RMSD 0.0295 Å | pass |
| `04_orbitals.py` rerun (Mes-Acr⁺, eosin Y) | localisation percentages identical to stored values | pass |
| MMFF vs DFT geometry, 6 molecules (`05_geometry_check.py`) | rigid: \|Δgap\| ≤ 0.03 eV; flexible/charged: \|Δgap\| ≤ 0.15 eV | documented limitation for the 10 MMFF-only molecules |
| Water / benzene optimiser check (commit `d07776c`) | O–H 0.967 Å, C–C 1.398 Å | pass (earlier work) |

## Not covered yet

These were planned in the previous version of this repo (commit `d07776c`) and are not part of the current results:

- TDDFT excitation energies and def2-TZVP single points (current orbitals are B3LYP/def2-SVP).
- D3(BJ) re-optimisation of the flexible dyes (eosin Y, riboflavin, rhodamine B); D3 is available, but those geometries are still without it.
- DFT geometries for the 10 MMFF-only molecules, most importantly RB²⁻ and FL²⁻ (eosin Y scaffold) and the large 4CzIPN / TPP.
- Anthraquinone-2-sulfonate, Ru(bpy)₃²⁺, fac-Ir(ppy)₃, decatungstate (TBADT) and GaN. The old
  `data/catalysts.csv` with these entries is still in git history.

## Log

- **2026-10-16** — Environment, catalyst list v1, representation step, optimiser validated on water/benzene (commit `d07776c`).
- **2026-10-17** — Full pipeline:
  - 16 organic photocatalysts: MMFF geometries; DFT optimisation of 3 molecules; B3LYP frontier
    orbitals for all 16; orbital localisation for Mes-Acr⁺ and eosin Y.
  - Solids (13 materials, Materials Project, anatase PBE bands) also done; not yet in this repo.
- **2026-10-18** — Organic part reorganised around scripts and verified by rerunning steps 1–4.
  Old `data/catalysts.csv` and `solids/qe/` kept until the solid-state update.
- **2026-10-24** — Review follow-up:
  - Removed the water redox lines from the orbital diagram and renamed the −HOMO/−LUMO−4.44 columns `*_qualitative`.
  - DFT-optimised eosin Y²⁻, riboflavin and rhodamine B⁺ (`--criteria flat`, `--start` added to step 2).
  - Added step 5, which quantifies the MMFF-vs-DFT geometry error. Eosin Y localisation recomputed on its DFT geometry.
