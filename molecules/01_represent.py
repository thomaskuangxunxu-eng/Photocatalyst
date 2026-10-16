"""Step 1 of the molecule pipeline: SMILES -> identifiers, 2D picture, 3D geometry, .xyz.

Reads data/catalysts.csv, keeps rows that have a SMILES, and for each molecule:
  * checks the SMILES parses and its formula matches `expected_formula` (catches typos)
  * writes canonical SMILES, InChI, InChIKey, formula, heavy-atom/total-atom counts
  * embeds 3D coordinates (ETKDGv3) and pre-optimizes with MMFF94 (UFF fallback)
  * writes molecules/xyz/<slug>.xyz (charge in the comment line)
Outputs: data/molecules_index.csv, results/figures/molecules_grid.png

Run from the repo root:  python molecules/01_represent.py
"""
from pathlib import Path
import re

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

ROOT = Path(__file__).resolve().parents[1]
XYZ_DIR = ROOT / "molecules" / "xyz"
FIG_DIR = ROOT / "results" / "figures"


def slugify(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def norm_formula(f: str) -> str:
    """RDKit writes charges as 'C20H6Br4O5-2' / 'C23H22N+'; normalise for comparison."""
    return f.replace("-1", "-").replace("+1", "+").strip()


def embed(mol: Chem.Mol, seed: int = 42) -> tuple[Chem.Mol, str]:
    molH = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(molH, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(molH, params) != 0:
            raise RuntimeError("3D embedding failed")
    if AllChem.MMFFHasAllMoleculeParams(molH):
        AllChem.MMFFOptimizeMolecule(molH, maxIters=2000)
        ff = "MMFF94"
    else:  # e.g. iodine-containing dyes
        AllChem.UFFOptimizeMolecule(molH, maxIters=2000)
        ff = "UFF"
    return molH, ff


def main() -> None:
    XYZ_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(ROOT / "data" / "catalysts.csv")
    df = df[df["smiles"].notna()].copy()

    rows, mols, legends = [], [], []
    for _, r in df.iterrows():
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            print(f"[FAIL] {r['name']}: SMILES did not parse")
            continue
        formula = CalcMolFormula(mol)
        charge = Chem.GetFormalCharge(mol)
        ok_formula = norm_formula(formula) == norm_formula(str(r["expected_formula"]))
        ok_charge = charge == int(r["charge"])
        molH, ff = embed(mol)
        slug = slugify(r["name"])
        xyz = Chem.MolToXYZBlock(molH).splitlines()
        xyz[1] = f"{r['name']} charge={charge} spin=0 geom={ff}"
        (XYZ_DIR / f"{slug}.xyz").write_text("\n".join(xyz) + "\n")

        rows.append({
            "name": r["name"], "slug": slug,
            "canonical_smiles": Chem.MolToSmiles(mol),
            "inchi": Chem.MolToInchi(mol),
            "inchikey": Chem.MolToInchiKey(mol),
            "formula": formula, "charge": charge,
            "n_heavy": mol.GetNumHeavyAtoms(), "n_atoms": molH.GetNumAtoms(),
            "n_electrons": sum(a.GetAtomicNum() for a in molH.GetAtoms()) - charge,
            "preopt_ff": ff, "formula_ok": ok_formula, "charge_ok": ok_charge,
        })
        flag = "OK  " if ok_formula and ok_charge else "FAIL"
        print(f"[{flag}] {r['name']:<28} {formula:<16} atoms={molH.GetNumAtoms():>3}  {ff}")
        mols.append(mol)
        legends.append(r["name"])

    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "data" / "molecules_index.csv", index=False)
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(320, 260), legends=legends)
    img.save(FIG_DIR / "molecules_grid.png")
    n_bad = (~(out.formula_ok & out.charge_ok)).sum()
    print(f"\n{len(out)} molecules written; {n_bad} failed the formula/charge check")


if __name__ == "__main__":
    main()
