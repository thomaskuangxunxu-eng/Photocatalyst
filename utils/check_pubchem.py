"""Cross-check every InChIKey in data/molecules_index.csv against PubChem.

Run locally (needs internet):  python utils/check_pubchem.py
For each molecule it prints the PubChem CID and IUPAC/title, so you can confirm the
SMILES encodes the compound you meant. Charged species (dye dianions, cations) often
have no exact PubChem entry; for those the script also tries the neutral-parent
skeleton (first 14 characters of the InChIKey).
"""
from pathlib import Path
import time

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def lookup(inchikey: str) -> str:
    r = requests.get(f"{BASE}/compound/inchikey/{inchikey}/property/Title,MolecularFormula/JSON", timeout=20)
    if r.ok:
        p = r.json()["PropertyTable"]["Properties"][0]
        return f"exact  CID {p['CID']}: {p.get('Title', '?')} ({p['MolecularFormula']})"
    r = requests.get(f"{BASE}/compound/inchikey/{inchikey[:14]}/property/Title/JSON", timeout=20)
    if r.ok:
        titles = [p.get("Title", "?") for p in r.json()["PropertyTable"]["Properties"][:3]]
        return "skeleton match: " + "; ".join(titles)
    return "NOT FOUND - check the SMILES"


def main() -> None:
    df = pd.read_csv(ROOT / "data" / "molecules_index.csv")
    for _, row in df.iterrows():
        print(f"{row['name']:<28} {lookup(row['inchikey'])}")
        time.sleep(0.25)  # PubChem asks for <= 5 requests/s


if __name__ == "__main__":
    main()
