"""Step 1: catalogue -> identifiers check, 2D grid, MMFF 3D geometry (.xyz), 3D viewer.

Reads data/organic_photocatalysts.csv. For each molecule:
  * parses the SMILES and checks formula + formal charge against the catalogue columns
  * recomputes the InChIKey and checks it matches
  * ETKDGv3 multi-conformer embedding (8 conformers, 30 if >3 rotatable bonds), MMFF94s,
    keeps the lowest-energy conformer -> molecules/xyz/<slug>.xyz (charge in comment line)
Outputs: molecules/xyz/*.xyz, results/molecules/figures/organic_photocatalysts_grid.png,
         results/molecules/organic_photocatalysts_3d.html (py3Dmol viewer, needs internet for 3Dmol.js)

Run from the repo root:  python molecules/01_represent.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, inchi, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

from utils.common import ROOT, write_xyz

cat = pd.read_csv(ROOT / "data/organic_photocatalysts.csv")
xyz_dir = ROOT / "molecules/xyz"; xyz_dir.mkdir(parents=True, exist_ok=True)
fig_dir = ROOT / "results/molecules/figures"; fig_dir.mkdir(parents=True, exist_ok=True)

rows, mols2d, legends, problems = [], [], [], []
for _, r in cat.iterrows():
    m = Chem.MolFromSmiles(r.smiles)
    if m is None:
        problems.append(f"{r.slug}: SMILES does not parse"); continue
    formula, q, key = rdMolDescriptors.CalcMolFormula(m), Chem.GetFormalCharge(m), inchi.MolToInchiKey(m)
    if formula != r.formula: problems.append(f"{r.slug}: formula {formula} != {r.formula}")
    if q != int(r.charge):   problems.append(f"{r.slug}: charge {q} != {r.charge}")
    if key != r.inchikey:    problems.append(f"{r.slug}: InChIKey {key} != {r.inchikey}")

    mh = Chem.AddHs(m)
    ps = AllChem.ETKDGv3(); ps.randomSeed = 0xF00D; ps.pruneRmsThresh = 0.5
    nconf = 30 if rdMolDescriptors.CalcNumRotatableBonds(mh) > 3 else 8
    cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=nconf, params=ps))
    res = AllChem.MMFFOptimizeMoleculeConfs(mh, mmffVariant="MMFF94s", maxIters=2000)
    emin, best = min((e, c) for (conv, e), c in zip(res, cids))
    syms = [a.GetSymbol() for a in mh.GetAtoms()]
    write_xyz(xyz_dir / f"{r.slug}.xyz", syms, mh.GetConformer(best).GetPositions(),
              f"{r['name']} | {formula} | charge={q} mult=1 | MMFF94s, lowest of {len(cids)} "
              f"ETKDGv3 confs | E_MMFF={emin:.2f} kcal/mol")
    rows.append(dict(slug=r.slug, n_conf=len(cids), all_mmff_converged=all(c == 0 for c, _ in res),
                     E_mmff_kcal=round(emin, 2), n_atoms=len(syms)))
    AllChem.Compute2DCoords(m); mols2d.append(m)
    qs = {1: "+1", -2: "-2", 0: "neutral"}.get(q, str(q))
    legends.append(f"{r['name'].split(' (')[0]}  ({r.abbrev})\n{formula}, q={qs}, lmax={r.lambda_max_nm} nm")

d = rdMolDraw2D.MolDraw2DCairo(1800, 1500, 450, 375)
d.drawOptions().legendFontSize = 20; d.drawOptions().addStereoAnnotation = False
d.DrawMolecules(mols2d, legends=legends); d.FinishDrawing()
(fig_dir / "organic_photocatalysts_grid.png").write_bytes(d.GetDrawingText())

# self-contained py3Dmol viewer (loads 3Dmol.js from the web)
xyzs = {s: (xyz_dir / f"{s}.xyz").read_text() for s in cat.slug}
meta = {r.slug: f"{r['name']} - {r.family}; {r.role}" for _, r in cat.iterrows()}
opts = "\n".join(f'<option value="{r.slug}">{r.abbrev} - {r["name"]}</option>' for _, r in cat.iterrows())
html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>body{{font-family:-apple-system,Helvetica,sans-serif;margin:18px}}
#v{{width:760px;height:560px;position:relative;border:1px solid #ccc}}
select{{font-size:15px;padding:4px}} #info{{margin-top:8px;color:#333;max-width:760px}}</style></head><body>
<h2>Organic photocatalysts - MMFF94s 3D structures</h2>
<select id="sel">{opts}</select> <label><input type="checkbox" id="stick"> stick only</label>
<div id="v"></div><div id="info"></div>
<script>
const XYZ={json.dumps(xyzs)}; const META={json.dumps(meta)};
let v=$3Dmol.createViewer(document.getElementById('v'),{{backgroundColor:'white'}});
function show(){{const s=document.getElementById('sel').value;v.clear();v.addModel(XYZ[s],'xyz');
 const st=document.getElementById('stick').checked;
 v.setStyle({{}},st?{{stick:{{radius:0.14}}}}:{{stick:{{radius:0.12}},sphere:{{scale:0.25}}}});
 v.zoomTo();v.render();document.getElementById('info').textContent=META[s];}}
document.getElementById('sel').onchange=show;document.getElementById('stick').onchange=show;show();
</script></body></html>"""
(ROOT / "results/molecules/organic_photocatalysts_3d.html").write_text(html)

print(pd.DataFrame(rows).to_string(index=False))
print("checks:", "all passed" if not problems else problems)
sys.exit(1 if problems else 0)
