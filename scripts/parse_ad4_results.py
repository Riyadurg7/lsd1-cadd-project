"""Parse AutoDock4 .dlg files into a results table."""
import os, re, glob, math
import pandas as pd

AD4  = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ad4")
DLG  = os.path.join(AD4, "dlg")
MINI = os.path.join(AD4, "..", "docking", "docking_set_minimal.csv")
VINA = os.path.join(AD4, "..", "docking", "docking_results_minimal.csv")

RT = 0.0019872 * 298.15          # kcal/mol at 298.15 K

mini = pd.read_csv(MINI)
meta = {r["lig_id"]: r for _, r in mini.iterrows()}

ROW = re.compile(r"^\s*(\d+)\s*\|\s*(-?\d+\.\d+)\s*\|\s*(\d+)\s*\|"
                 r"\s*(-?\d+\.\d+)\s*\|\s*(\d+)\s*\|", re.M)

rows = []
for path in sorted(glob.glob(os.path.join(DLG, "*.dlg"))):
    lid = os.path.splitext(os.path.basename(path))[0]
    txt = open(path, errors="replace").read()

    start = txt.find("CLUSTERING HISTOGRAM")
    end   = txt.find("RMSD TABLE", start)
    block = txt[start:end] if start >= 0 and end > start else ""

    clusters = [(int(m.group(1)), float(m.group(2)), int(m.group(5)))
                for m in ROW.finditer(block)]

    best_e = clusters[0][1] if clusters else None
    top_pop = clusters[0][2] if clusters else None
    n_clust = len(clusters) if clusters else None

    # how many of the 10 runs landed in the best-energy cluster
    ki_nM = round(math.exp(best_e / RT) * 1e9, 1) if best_e is not None else None

    r = meta.get(lid, {})
    rows.append({
        "lig_id": lid,
        "group": r.get("group", ""),
        "AD4_energy": best_e,
        "Ki_nM": ki_nM,
        "n_clusters": n_clust,
        "runs_in_best": top_pop,
        "Tanimoto": r.get("Tanimoto", ""),
        "pIC50": r.get("pIC50", ""),
    })

df = pd.DataFrame(rows)

vina = pd.read_csv(VINA)[["lig_id", "score"]].rename(columns={"score": "Vina"})
df = df.merge(vina, on="lig_id", how="left")
df = df.sort_values("AD4_energy").reset_index(drop=True)
df.index += 1

out = os.path.join(AD4, "autodock4_results.csv")
df.to_csv(out, index=False)

pd.set_option("display.width", 220)
print("AUTODOCK4 RESULTS - ranked by binding energy (most negative = best)")
print()
print(df[["lig_id", "group", "AD4_energy", "Ki_nM", "n_clusters",
          "runs_in_best", "Vina", "Tanimoto"]].to_string())

print("\nBY GROUP")
print(df.groupby("group")["AD4_energy"].agg(["count", "mean", "min"]).round(2).to_string())

both = df.dropna(subset=["AD4_energy", "Vina"])
if len(both) > 2:
    from scipy.stats import pearsonr, spearmanr
    r1, p1 = pearsonr(both["AD4_energy"], both["Vina"])
    r2, p2 = spearmanr(both["AD4_energy"], both["Vina"])
    print(f"\nAUTODOCK4 vs VINA   Pearson r = {r1:.3f} (p = {p1:.4f})"
          f"   Spearman rho = {r2:.3f} (p = {p2:.4f})")

print("\nsaved:", out)
