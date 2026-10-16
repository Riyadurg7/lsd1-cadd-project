"""Re-prepare ligands that contain Meeko macrocycle glue atoms (G0/CG0),
which AutoDock4 cannot read. Rigid macrocycles are used instead."""
import os, sys, inspect
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # project root
LIGD = os.path.join(HERE, "docking", "ligands")

TARGETS = sys.argv[1:] or ["E0008"]

sel = pd.read_csv(os.path.join(HERE, "enamine", "enamine_selected.csv"))
mini = pd.read_csv(os.path.join(HERE, "docking", "docking_set_minimal.csv"))
smiles_by_id = dict(zip(mini["lig_id"], mini["clean_smiles"]))

# find the keyword this Meeko build uses for rigid macrocycles
sig = inspect.signature(MoleculePreparation.__init__).parameters
if "rigid_macrocycles" in sig:
    kw = {"rigid_macrocycles": True}
elif "macrocycle" in sig:
    kw = {"macrocycle": False}
else:
    kw = {}
print("MoleculePreparation kwargs:", kw)

for lid in TARGETS:
    smi = smiles_by_id.get(lid)
    if smi is None:
        print(lid, "- not in the minimal set, skipped")
        continue

    m = Chem.MolFromSmiles(smi)
    frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        m = max(frags, key=lambda x: x.GetNumHeavyAtoms())
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=400)
    except Exception:
        pass

    setups = MoleculePreparation(**kw).prepare(m)
    txt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        print(lid, "- write failed:", err)
        continue

    path = os.path.join(LIGD, lid + ".pdbqt")
    open(path, "w").write(txt)

    types = set()
    for line in txt.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            types.add(line.split()[-1])
    glue = {t for t in types if t.startswith(("G", "CG")) and t not in ("C",)}
    print(f"{lid}: rewritten | types = {sorted(types)} | glue atoms = {sorted(glue) or 'none'}")
