# Photocatalyst catalogue and electronic-structure workflow — results report

_Generated 2026-10-17; geometry section and redox caveat updated 2026-10-24. Paths are relative to the repository root. How to rerun each step is in [README.md](README.md)._

## 1. Organic photocatalysts (16)

Catalogue: [data/organic_photocatalysts.csv](data/organic_photocatalysts.csv). Structure grid: [results/molecules/figures/organic_photocatalysts_grid.png](results/molecules/figures/organic_photocatalysts_grid.png).

**Representation.** Each entry carries a canonical RDKit SMILES (concise, human-editable, encodes charge), a standard InChIKey
(fixed-length hash for database lookup/deduplication), and a PubChem CID. Names were resolved against PubChem. For 9 of 16 entries the stored SMILES/InChIKey matches the PubChem record exactly. For the other 7, the column `pubchem_parent_cid` points to a parent record, not the stored species: for salts the counter-ion was stripped and the active ion kept (EY²⁻, RB²⁻, RhB⁺, MB⁺, Mes-Acr⁺); FL²⁻ is a quinoid dianion built locally, while CID 16850 is the neutral lactone; Me-PDI has no PubChem record and was built from PTCDA (CID 67191). See `pubchem_cid_relation`.

**3D structures / XYZ.** RDKit ETKDGv3 multi-conformer embedding (8 conformers, 30 for >3 rotatable bonds, RMS prune 0.5 Å),
MMFF94s optimisation, lowest-energy conformer written to XYZ with charge/multiplicity in the comment line
([molecules/xyz/](molecules/xyz/)). Interactive py3Dmol viewer with a molecule selector:
[results/molecules/organic_photocatalysts_3d.html](results/molecules/organic_photocatalysts_3d.html). Static montage: [results/molecules/figures/structures_3d_montage.png](results/molecules/figures/structures_3d_montage.png).

**DFT geometry optimisation** (PySCF 2.14.0 + geomeTRIC), B3LYP/def2-SVP, density fitting, ddCOSMO water (ε = 78.36):
[results/molecules/geomopt_summary.csv](results/molecules/geomopt_summary.csv), geometries [molecules/xyz_opt/](molecules/xyz_opt/),
overlay [results/molecules/figures/geomopt_overlay_benzophenone.png](results/molecules/figures/geomopt_overlay_benzophenone.png).
Six molecules were optimised: the three rigid aryl ketones/quinones, plus the flexible or charged eosin Y²⁻, riboflavin and rhodamine B⁺.

| molecule | basis fns | steps | criteria | ΔE vs MMFF geom. (kcal/mol) | heavy-atom RMSD (Å) | max atom dev. (Å) |
|---|---|---|---|---|---|---|
| benzophenone | 246 | 11 | gau | -1.27 | 0.0956 | 0.2625 |
| anthraquinone | 264 | 15 | gau | -1.42 | 0.0192 | 0.0367 |
| thioxanthone | 254 | 10 | gau | -1.88 | 0.0295 | 0.0618 |
| eosinY | 508 | 31 | flat* | -17.73 | 0.3163 | 1.2972 |
| riboflavin | 478 | 18 | flat* | -11.13 | 0.1732 | 0.6507 |
| rhodamineB | 617 | 20 | flat* | -20.23 | 0.3983 | 1.2664 |

\* `flat` criteria: same gradient thresholds as `gau` (RMS/max 3×10⁻⁴/4.5×10⁻⁴ Ha/bohr), ΔE < 10⁻⁵ Ha, no displacement thresholds.
Eosin Y²⁻ met the gradient thresholds after about 20 `gau` steps, but a very flat benzoate/carboxylate torsion kept the displacement test
from settling. The run was resumed from step 27 under `flat` criteria. Riboflavin and rhodamine B⁺ (ribityl chain, NEt₂ groups) were run
under `flat` from the start. Rhodamine B⁺ took 20 steps and ≈19 h wall time on this machine.

**How much do the MMFF geometries matter?** [results/molecules/geometry_sensitivity.csv](results/molecules/geometry_sensitivity.csv)
(`molecules/05_geometry_check.py`). Single points on the MMFF geometry versus the DFT geometry, Δ = DFT − MMFF:

| PC | q | heavy-atom RMSD (Å) | ΔHOMO (eV) | ΔLUMO (eV) | Δgap (eV) |
|---|---|---|---|---|---|
| AQ | 0 | 0.0192 | -0.026 | -0.02 | 0.006 |
| TX | 0 | 0.0295 | 0.024 | -0.002 | -0.027 |
| BP | 0 | 0.0956 | -0.089 | -0.087 | 0.001 |
| RF | 0 | 0.1732 | 0.066 | 0.014 | -0.052 |
| EY2- | -2 | 0.3163 | -0.005 | 0.143 | 0.149 |
| RhB+ | 1 | 0.3983 | -0.069 | -0.003 | 0.065 |

- **Rigid molecules** (BP, AQ, TX): the gap changes by ≤ 0.03 eV. The absolute orbital energies shift by up to
  0.09 eV (benzophenone, whose phenyl rings untwist by about 6°).
- **Flexible / charged molecules:** the geometries change substantially (heavy-atom RMSD 0.17–0.40 Å,
  -11 to -20 kcal/mol). The gap shifts by up to +0.15 eV (EY2-), mainly through its LUMO.
- **Size of the error:** a geometry error of about 0.15 eV is a third of the 0.45 eV scatter against experiment.
  It is small enough not to change any qualitative conclusion here, but too large to ignore for quantitative work.
- **Still on MMFF geometries (†, 10 of 16):** DCA, BZ, PTH, RB2-, FL2-, MB+, Mes-Acr+, Me-PDI, TPP, 4CzIPN. RB²⁻ and FL²⁻ share eosin Y's benzoate-xanthene scaffold and
  probably carry a similar ~0.1–0.15 eV gap error. 4CzIPN, TPP, Mes-Acr⁺ and MB⁺ are untested. Run `molecules/02_optimize.py <slug> --criteria flat`,
  then steps 3 and 5, to bring any of them onto a DFT geometry.

**Frontier orbitals** — B3LYP/def2-SVP (def2 ECP on I), DF, ddCOSMO water, all 16 SCFs converged
([results/molecules/frontier_orbitals.csv](results/molecules/frontier_orbitals.csv), diagram [results/molecules/figures/homo_lumo_diagram.png](results/molecules/figures/homo_lumo_diagram.png)).
DFT geometries for TX, BP, AQ, EY2-, RF, RhB+; MMFF geometries (†) for the rest.

| PC | family | q | geometry | HOMO (eV) | LUMO (eV) | gap (eV) | lowest abs. band (nm) | −LUMO−4.44 (qualitative index, not V vs SHE) |
|---|---|---|---|---|---|---|---|---|
| TPP | porphyrin | 0 | MMFF † | -5.1229 | -3.4268 | 1.6961 | 647 | -1.01 |
| Me-PDI | rylene diimide | 0 | MMFF † | -5.8143 | -3.5267 | 2.2876 | 525 | -0.91 |
| MB+ | phenothiazinium dye | 1 | MMFF † | -5.8746 | -3.4453 | 2.4293 | 664 | -0.99 |
| RB2- | xanthene dye (heavy-atom) | -2 | MMFF † | -5.1405 | -2.4285 | 2.712 | 549 | -2.01 |
| DCA | cyanoarene | 0 | MMFF † | -6.0727 | -3.3068 | 2.7659 | 420 | -1.13 |
| RhB+ | rhodamine (xanthene) dye | 1 | DFT | -5.7773 | -2.9474 | 2.8299 | 554 | -1.49 |
| FL2- | xanthene dye | -2 | MMFF † | -4.6206 | -1.7525 | 2.8681 | 490 | -2.69 |
| EY2- | xanthene dye | -2 | DFT | -5.029 | -2.0747 | 2.9544 | 539 | -2.37 |
| 4CzIPN | cyanoarene / TADF donor-acceptor | 0 | MMFF † | -5.85 | -2.7794 | 3.0706 | 435 | -1.66 |
| Mes-Acr+ | acridinium | 1 | MMFF † | -6.7851 | -3.6183 | 3.1668 | 430 | -0.82 |
| RF | flavin | 0 | DFT | -6.3163 | -2.9415 | 3.3748 | 445 | -1.5 |
| BZ | 1,2-diketone | 0 | MMFF † | -6.9214 | -3.0754 | 3.846 | 370 | -1.36 |
| TX | aryl ketone (thio) | 0 | DFT | -6.1212 | -2.1486 | 3.9726 | 380 | -2.29 |
| PTH | phenothiazine | 0 | MMFF † | -4.9883 | -0.8909 | 4.0974 | 316 | -3.55 |
| AQ | quinone | 0 | DFT | -7.3138 | -3.0409 | 4.2729 | 325 | -1.4 |
| BP | aryl ketone | 0 | DFT | -6.9985 | -2.0343 | 4.9643 | 345 | -2.41 |

The Kohn–Sham gap correlates with the lowest experimental absorption band (r = 0.89, MAE 0.46 eV, n = 16), but it is not an
excitation energy. Largest outlier: benzophenone (+1.37 eV), whose weak first band is n→π*. TPP is compared against its Q_x(0,0) band
(647 nm), not the Soret (419 nm, S₂). For excitation energies use TD-DFT.

**The orbital energies are not redox potentials.** The last column (−LUMO − 4.44) is only a rough ordering index. It ignores the solvation
free energy of the ion and the relaxation on reduction, and for the dianions it lands more than 1 V away from measured potentials
(e.g. EY²⁻: -2.37 here vs ≈ −0.8 V vs SHE for EY/EY•⁻).
For this reason the diagram does not draw the H⁺/H₂ or O₂/H₂O levels, and the CSV columns are named `*_qualitative`.

**Orbital localisation** ([results/molecules/figures/frontier_orbitals_isosurfaces.png](results/molecules/figures/frontier_orbitals_isosurfaces.png), cubes [results/molecules/cubes/](results/molecules/cubes/),
table [results/molecules/orbital_localization.csv](results/molecules/orbital_localization.csv)). Mulliken orbital populations: Mes-Acr⁺ LUMO is 98 % on acridinium; the
mesityl donor pair (HOMO−1/HOMO−2, 96–100 % mesityl) is within 13 meV of an acridinium-centred HOMO, so the KS ordering among
these three is not meaningful at this level. Eosin Y²⁻ (DFT geometry) HOMO and LUMO are 98 % and 96 % on the xanthene core (π→π*); on the MMFF geometry the values were 91 % and 97 %.

## 2. Tool roles and limits

- RDKit: identifiers, 2D depiction, conformers, XYZ export. MMFF94s geometries are adequate for rigid chromophores (gap within 0.03 eV of
  the DFT-geometry value). For flexible or charged dyes they shift the gap by up to 0.15 eV, so optimise those with DFT.
- PySCF: molecular DFT (single points, geometry optimisation with geomeTRIC, cube files).
- py3Dmol: interactive inspection in HTML. It has no headless PNG export, so static isosurfaces were rendered by marching cubes in matplotlib.

**Dispersion.** All stored results are without dispersion correction. The PyPI `pyscf-dispersion` 1.0.0 cannot load on Apple
Silicon: it calls `load_library()` without its path argument and ships x86_64-only binaries. D3/D4 now run through
`utils/dispersion_shim.py`, which sits on the conda-forge `dftd3-python` / `dftd4-python` builds of the same libraries.
D3 does not change the Kohn–Sham orbitals at a fixed geometry. It affects the results only through optimised geometries:
re-optimising benzophenone and thioxanthone with D3(BJ) moved them by ≤ 0.022 Å heavy-atom RMSD and their gaps by ≤ 0.006 eV.
The effect on the flexible dyes is not yet measured. Ground-state KS orbital energies are not redox potentials (see the frontier-orbital section).
