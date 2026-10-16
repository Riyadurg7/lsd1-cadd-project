"""Measure how close each docked pose comes to the FAD cofactor."""
import glob, os
import numpy as np

AD4 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ad4")

rec = open(os.path.join(AD4, "receptor.pdb")).read().splitlines()
fad = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                for l in rec
                if l.startswith("HETATM") and l[17:20].strip() == "FAD"])
print("FAD atoms in receptor:", len(fad))
print()
print("{:<10}{:>20}{:>20}".format("ligand", "min dist to FAD", "atoms < 4.5 A"))
print("-" * 50)

for p in sorted(glob.glob(os.path.join(AD4, "complexes", "*_pose.pdb"))):
    lid = os.path.basename(p).split("_")[0]
    lig = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    for l in open(p) if l.startswith("HETATM")])
    d = np.sqrt(((lig[:, None, :] - fad[None, :, :]) ** 2).sum(-1))
    print("{:<10}{:>20.2f}{:>20d}".format(lid, d.min(), int((d < 4.5).sum())))

print()
print("Interpretation:")
print("  < 4.0 A  direct contact with the cofactor")
print("  4-6 A    adjacent, in the catalytic cleft")
print("  > 8 A    binding elsewhere in the substrate channel")
