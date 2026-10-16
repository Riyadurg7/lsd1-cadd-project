"""
Protein-ligand interaction profiling (PLIP) for the top docked poses.

Note: the pip Open Babel wheel is built without the InChI plugin, so PLIP's
call to write(format='inchikey') raises. The InChIKey is only an identifier
in the report - it plays no part in interaction detection - so it is patched
to return an empty string.
"""
import os, glob
import pandas as pd

# ---- patch out the InChIKey call before PLIP uses it ------------------
import openbabel.pybel as pybel

_orig_write = pybel.Molecule.write


def _safe_write(self, format="smi", *a, **kw):
    if format == "inchikey":
        return ""
    return _orig_write(self, format, *a, **kw)


pybel.Molecule.write = _safe_write

from plip.structure.preparation import PDBComplex   # noqa: E402

AD4 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ad4")
CPX = os.path.join(AD4, "complexes")

res = pd.read_csv(os.path.join(AD4, "autodock4_results.csv"))
energy = dict(zip(res["lig_id"], res["AD4_energy"]))
group = dict(zip(res["lig_id"], res["group"]))

rows = []
for path in sorted(glob.glob(os.path.join(CPX, "complex_*.pdb"))):
    lid = os.path.basename(path)[8:-4]
    try:
        cx = PDBComplex()
        cx.load_pdb(path)
        cx.analyze()
    except Exception as e:
        print(f"{lid}: PLIP failed - {type(e).__name__}: {e}")
        continue

    before = len(rows)

    def collect(kind, objs, dattr=None):
        out = []
        for o in objs:
            d = getattr(o, dattr, None) if dattr else None
            out.append({
                "ligand": lid,
                "group": group.get(lid, ""),
                "energy": energy.get(lid),
                "interaction": kind,
                "residue": f"{getattr(o, 'restype', '?')}{getattr(o, 'resnr', '')}",
                "distance": round(d, 2) if isinstance(d, float) else "",
            })
        return out

    for key, site in cx.interaction_sets.items():
        if "LIG" not in key:
            continue
        rows += collect("H-bond", list(site.hbonds_pdon) + list(site.hbonds_ldon),
                        "distance_ad")
        rows += collect("Hydrophobic", site.hydrophobic_contacts, "distance")
        rows += collect("pi-stacking", site.pistacking, "distance")
        rows += collect("Salt bridge",
                        list(site.saltbridge_lneg) + list(site.saltbridge_pneg),
                        "distance")
        rows += collect("Water bridge", site.water_bridges, "distance_aw")
        rows += collect("Halogen bond", site.halogen_bonds, "distance")

    print(f"{lid}: {len(rows) - before} interactions")

df = pd.DataFrame(rows)
if df.empty:
    print("\nNo interactions parsed.")
    raise SystemExit

out = os.path.join(AD4, "plip_interactions.csv")
df.to_csv(out, index=False)

pd.set_option("display.width", 200)
print("\n" + "=" * 70)
print("INTERACTION COUNTS PER LIGAND")
print("=" * 70)
print(pd.crosstab(df["ligand"], df["interaction"]).to_string())

print("\n" + "=" * 70)
print("CONTACTS WITH THE FAD COFACTOR")
print("=" * 70)
fad = df[df["residue"].str.contains("FAD", na=False)]
print(fad.to_string(index=False) if len(fad) else "none detected")

print("\n" + "=" * 70)
print("MOST FREQUENTLY CONTACTED RESIDUES")
print("=" * 70)
print(df["residue"].value_counts().head(12).to_string())

print("\n" + "=" * 70)
print("PER-LIGAND DETAIL")
print("=" * 70)
for lid in sorted(df["ligand"].unique(), key=lambda x: energy.get(x, 0)):
    sub = df[df["ligand"] == lid]
    print(f"\n{lid}   {group.get(lid,'')}   {energy.get(lid)} kcal/mol")
    for kind in sub["interaction"].unique():
        resl = sorted(set(sub[sub["interaction"] == kind]["residue"]))
        print(f"   {kind:<14} {', '.join(resl)}")

print("\nsaved:", out)
