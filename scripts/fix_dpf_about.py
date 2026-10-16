"""
Fix the 'about' line in every AutoDock4 dpf.

'about' is the ligand's own centre of coordinates in its PDBQT file - the point
AutoDock rotates around. It is NOT the grid centre. Setting it to the grid
centre translates the ligand far outside the box, which AutoDock reports as a
huge positive energy.
"""
import os, glob

AD4 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ad4")


def centroid(pdbqt):
    xs = ys = zs = 0.0
    n = 0
    for line in open(pdbqt):
        if line.startswith(("ATOM", "HETATM")):
            xs += float(line[30:38])
            ys += float(line[38:46])
            zs += float(line[46:54])
            n += 1
    return (xs / n, ys / n, zs / n) if n else (0.0, 0.0, 0.0)


fixed = 0
for dpf in sorted(glob.glob(os.path.join(AD4, "*.dpf"))):
    lid = os.path.splitext(os.path.basename(dpf))[0]
    lig = os.path.join(AD4, lid + ".pdbqt")
    if not os.path.exists(lig):
        continue

    cx, cy, cz = centroid(lig)
    out = []
    for line in open(dpf):
        if line.startswith("about "):
            out.append(f"about {cx:.3f} {cy:.3f} {cz:.3f}\n")
        else:
            out.append(line)
    open(dpf, "w", newline="\n").writelines(out)
    fixed += 1
    if fixed <= 5:
        print(f"{lid}: about {cx:.3f} {cy:.3f} {cz:.3f}")

print(f"\n{fixed} dpf files corrected")
