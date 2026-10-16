"""
Select Enamine similarity hits and prepare them for docking.

  - keeps every hit at or above STRICT (default 0.90)
  - clusters the band between LOOSE and STRICT, taking one representative each
  - writes 3D PDBQT files into docking/ligands/
  - appends the selected compounds to docking/docking_set_prepared.csv
"""
import argparse, os, sys
import pandas as pd
import numpy as np

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # project root
HITS = os.path.join(HERE, "enamine", "Enamine_similarity_hits.csv")
DOCK = os.path.join(HERE, "docking")
LIGD = os.path.join(DOCK, "ligands")
SETF = os.path.join(DOCK, "docking_set_prepared.csv")

ap = argparse.ArgumentParser()
ap.add_argument("--strict", type=float, default=0.90)
ap.add_argument("--loose",  type=float, default=0.85)
ap.add_argument("--n-tier2", type=int, default=50,
                help="how many clustered representatives from the loose band")
ap.add_argument("--cluster-cutoff", type=float, default=0.40,
                help="Butina distance cutoff (larger = fewer, more diverse clusters)")
args = ap.parse_args()

for p in (HITS, SETF, LIGD):
    if not os.path.exists(p):
        sys.exit("ERROR: missing " + p)

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

h = pd.read_csv(HITS).drop_duplicates("Enamine_SMILES")
print(f"hits loaded: {len(h):,}")

# ---- tier 1: everything at or above the strict cutoff -----------------
tier1 = h[h["Tanimoto"] >= args.strict].sort_values("Tanimoto", ascending=False)
print(f"tier 1  (>= {args.strict}): {len(tier1)}")

# ---- tier 2: cluster the band, take one representative per cluster ----
band = h[(h["Tanimoto"] >= args.loose) & (h["Tanimoto"] < args.strict)]
band = band.sort_values("Tanimoto", ascending=False).head(400).reset_index(drop=True)
print(f"band    ({args.loose}-{args.strict}): {len(band)} considered")

fps, keep_idx = [], []
for i, smi in enumerate(band["Enamine_SMILES"]):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        continue
    fps.append(gen.GetFingerprint(m))
    keep_idx.append(i)
band = band.iloc[keep_idx].reset_index(drop=True)

dists = []
for i in range(1, len(fps)):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
    dists.extend(1.0 - np.array(sims))

clusters = Butina.ClusterData(dists, len(fps), args.cluster_cutoff, isDistData=True)
print(f"clusters in band: {len(clusters)}")

reps = [sorted(c, key=lambda j: -band.loc[j, "Tanimoto"])[0] for c in clusters]
reps = sorted(reps, key=lambda j: -band.loc[j, "Tanimoto"])[:args.n_tier2]
tier2 = band.iloc[reps]
print(f"tier 2  representatives: {len(tier2)}")

sel = pd.concat([tier1, tier2], ignore_index=True).drop_duplicates("Enamine_SMILES")
print(f"\nselected for docking: {len(sel)}")

# ---- prepare PDBQT ----------------------------------------------------
def make_pdbqt(smi, path):
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return False, "unparseable"
        frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
        if len(frags) > 1:
            m = max(frags, key=lambda x: x.GetNumHeavyAtoms())
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE) != 0:
            return False, "embed failed"
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=400)
        except Exception:
            pass
        setups = MoleculePreparation().prepare(m)
        if not setups:
            return False, "meeko empty"
        txt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return False, "write failed"
        open(path, "w").write(txt)
        return True, "ok"
    except Exception as e:
        return False, type(e).__name__

rows, failed = [], 0
for n, (_, r) in enumerate(sel.iterrows(), 1):
    lid = f"E{n:04d}"
    ok, why = make_pdbqt(r["Enamine_SMILES"], os.path.join(LIGD, lid + ".pdbqt"))
    if ok:
        rows.append({"clean_smiles": r["Enamine_SMILES"],
                     "pIC50": "",
                     "group": "enamine_hit",
                     "lig_id": lid,
                     "Tanimoto": r["Tanimoto"]})
    else:
        failed += 1

print(f"prepared {len(rows)} | failed {failed}")

new = pd.DataFrame(rows)
new[["Enamine_SMILES_dup"]] = new[["clean_smiles"]]
new.drop(columns=["Enamine_SMILES_dup"], inplace=True)
new.to_csv(os.path.join(HERE, "enamine", "enamine_selected.csv"), index=False)

# ---- append to the docking set ----------------------------------------
old = pd.read_csv(SETF)
old = old[old["group"] != "enamine_hit"]           # avoid duplicates on re-run
cols = [c for c in old.columns]
merged = pd.concat([old, new[[c for c in cols if c in new.columns]]],
                   ignore_index=True)
merged.to_csv(SETF, index=False)

print("\nDOCKING SET NOW:")
print(merged["group"].value_counts().to_string())
print("total:", len(merged))
print("\nwritten:", SETF)
