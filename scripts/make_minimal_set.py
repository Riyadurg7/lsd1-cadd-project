"""Build the minimal docking set: the 14 Enamine hits at >=0.90 plus 6 reference actives."""
import os
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # project root
DOCK = os.path.join(HERE, "docking")

sel = pd.read_csv(os.path.join(HERE, "enamine", "enamine_selected.csv"))
prep = pd.read_csv(os.path.join(DOCK, "docking_set_prepared.csv"))

# 14 Enamine hits at or above 0.90
ena = sel[sel["Tanimoto"] >= 0.90].copy()
ena["group"] = "enamine_hit"

# 6 most potent known actives, as a reference point
act = (prep[prep["group"] == "active"]
       .sort_values("pIC50", ascending=False)
       .head(6)
       .copy())
act["Tanimoto"] = ""

cols = ["lig_id", "group", "pIC50", "clean_smiles", "Tanimoto"]
mini = pd.concat([ena[cols], act[cols]], ignore_index=True)

# keep only ligands whose PDBQT actually exists
mini = mini[[os.path.exists(os.path.join(DOCK, "ligands", l + ".pdbqt"))
             for l in mini["lig_id"]]].reset_index(drop=True)

out = os.path.join(DOCK, "docking_set_minimal.csv")
mini.to_csv(out, index=False)

print("MINIMAL DOCKING SET")
print(mini["group"].value_counts().to_string())
print("total:", len(mini))
print("\nwritten:", out)
print()
print(mini[["lig_id", "group", "Tanimoto", "pIC50"]].to_string(index=False))
