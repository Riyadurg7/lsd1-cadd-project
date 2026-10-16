"""
Extract the lowest-energy docked pose from each AutoDock4 .dlg and build
receptor+ligand complex PDB files for PLIP.
"""
import os, re, glob
import pandas as pd

AD4 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ad4")
DLG = os.path.join(AD4, "dlg")
OUT = os.path.join(AD4, "complexes")
os.makedirs(OUT, exist_ok=True)

RES = pd.read_csv(os.path.join(AD4, "autodock4_results.csv"))
TOP = RES.nsmallest(5, "AD4_energy")["lig_id"].tolist()
print("top 5 by AutoDock4 energy:", TOP)

receptor = [l for l in open(os.path.join(AD4, "receptor.pdb")).read().splitlines()
            if l.startswith(("ATOM", "HETATM"))]
print("receptor atoms:", len(receptor))


def best_pose(dlg_path):
    """Return the ATOM/HETATM lines of the lowest-energy DOCKED conformation."""
    text = open(dlg_path, errors="replace").read()
    # each docked conformation is a MODEL block with DOCKED: prefixes
    blocks = re.findall(r"DOCKED: MODEL(.*?)DOCKED: ENDMDL", text, re.S)
    best, best_e = None, None
    for b in blocks:
        m = re.search(r"Estimated Free Energy of Binding\s*=\s*(-?\d+\.\d+)", b)
        if not m:
            continue
        e = float(m.group(1))
        if best_e is None or e < best_e:
            best_e, best = e, b
    if best is None:
        return None, None

    atoms = []
    for line in best.splitlines():
        line = line.replace("DOCKED: ", "", 1)
        if line.startswith(("ATOM", "HETATM")):
            atoms.append(line)
    return atoms, best_e


for lid in TOP:
    dlg = os.path.join(DLG, lid + ".dlg")
    atoms, energy = best_pose(dlg)
    if not atoms:
        print(f"{lid}: no pose found")
        continue

    # rewrite as HETATM, residue LIG, chain Z, and drop the PDBQT charge/type tail
    lig = []
    for i, l in enumerate(atoms, 1):
        x, y, z = l[30:38], l[38:46], l[46:54]
        elem = l[77:79].strip() if len(l) > 77 else l[12:14].strip()
        elem = re.sub(r"[^A-Za-z]", "", elem)[:2] or "C"
        if elem.upper() in ("HD", "HS"):
            elem = "H"
        elif elem.upper() in ("NA", "NS"):
            elem = "N"
        elif elem.upper() in ("OA", "OS"):
            elem = "O"
        elif elem.upper() == "SA":
            elem = "S"
        elif elem.upper() == "A":
            elem = "C"
        name = (elem + str(i))[:4]
        # strict PDB columns:
        #  1-6 record | 7-11 serial | 13-16 name | 17 altLoc | 18-20 resName
        #  22 chain | 23-26 resSeq | 31-38 x | 39-46 y | 47-54 z
        #  55-60 occ | 61-66 bfac | 77-78 element
        lig.append(
            "HETATM"
            + f"{i:5d}"
            + " "
            + f"{name:<4}"
            + " "
            + "LIG"
            + " "
            + "Z"
            + f"{901:4d}"
            + " " * 4
            + f"{x}{y}{z}"
            + "  1.00  0.00"
            + " " * 10
            + f"{elem.upper():>2}"
        )

    # ligand-only PDB
    open(os.path.join(OUT, f"{lid}_pose.pdb"), "w").write("\n".join(lig) + "\nEND\n")

    # complex
    cpx = os.path.join(OUT, f"complex_{lid}.pdb")
    open(cpx, "w").write("\n".join(receptor + lig) + "\nEND\n")

    print(f"{lid}: {len(lig)} ligand atoms, best pose {energy} kcal/mol -> {os.path.basename(cpx)}")

print("\ncomplexes written to:", OUT)
