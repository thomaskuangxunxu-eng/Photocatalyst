# photocatalyst-dft

Computational survey of reported photocatalysts: frontier orbitals (HOMO/LUMO, TDDFT) for organic
dyes with PySCF, then band gaps and band edges for solid photocatalysts via Materials Project and
periodic DFT (PySCF PBC, Quantum ESPRESSO).

## Setup

```bash
conda env create -f environment.yml
conda activate photocat
```

## Pipeline (molecules)

| Step | Command | Output |
| --- | --- | --- |
| 1. Represent | `python molecules/01_represent.py` | `data/molecules_index.csv` (canonical SMILES, InChI, InChIKey), `molecules/xyz/*.xyz`, `results/figures/molecules_grid.png` |
| 1b. Verify IDs | `python utils/check_pubchem.py` | PubChem title per InChIKey (needs internet) |
| 2. Optimize | `python molecules/02_optimize.py molecules/xyz/<name>.xyz [--solvent water]` | `results/opt/<name>/` (optimized xyz, `summary.json`) |
| 3. Frontier orbitals + TDDFT | `molecules/03_frontier.py` (Wed 21 Oct) | `results/molecules.csv` |

Default level of theory: B3LYP-D3(BJ)/def2-SVP for geometries; def2-TZVP single points for orbital
energies. Iodine and transition metals get def2 effective core potentials automatically.

## Data

`data/catalysts.csv`: 29 photocatalysts (14 organic with SMILES, 3 metal complexes/POMs for later,
12 inorganic solids). Each SMILES is the species that absorbs under typical conditions (e.g. the
eosin Y and fluorescein dianions), and its formula and charge are checked automatically against
`expected_formula`. Experimental λmax and band gaps are approximate literature values still to be
pinned to specific references (`reference` column).

## Sanity checks

| Check | Result | Reference | Status |
| --- | --- | --- | --- |
| All 14 SMILES parse; formula + charge match | 14/14 | `expected_formula` column | pass |
| Water, B3LYP-D3/def2-SVP opt | O–H 0.967 Å, ∠HOH 103.1° | exp. 0.958 Å, 104.5° (typical DZ-basis deviation) | pass |
| Benzene, B3LYP-D3/def2-SVP opt | C–C 1.398 Å, C–H 1.093 Å; 65 s on 2 cores | exp. C–C 1.397 Å, C–H 1.084 Å | pass |

## Log

- **2026-10-16** — Environment set up (PySCF 2.14, RDKit, geomeTRIC, pymatgen, mp-api). Catalyst list v1.
  Representation step done for 14 organic photocatalysts. Optimizer validated on water and benzene.
