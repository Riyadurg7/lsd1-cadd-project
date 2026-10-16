"""
Set up an AutoDock4 run for LSD1 using the course workflow
(autogrid4 -> affinity maps, autodock4 -> genetic-algorithm docking).

Creates D:\\DES\\CADD\\Project\\ad4\\ containing the receptor, the ligands,
a grid parameter file, one docking parameter file per ligand, and a bash
script that runs everything inside WSL.
"""
import os, shutil
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # project root
DOCK = os.path.join(HERE, "docking")
SRC  = os.path.join(HERE, "dockingflow")
AD4  = os.path.join(HERE, "ad4")

os.makedirs(AD4, exist_ok=True)
os.makedirs(os.path.join(AD4, "dlg"), exist_ok=True)

# ---- box, matching the Vina run -------------------------------------
CENTER  = (-2.722, 43.384, 85.953)
SPACING = 0.375
NPTS    = (60, 60, 60)          # 60 * 0.375 = 22.5 A per side

RECEPTOR_TYPES = "A C HD N NA OA P S"
LIGAND_TYPES   = "A C Cl F HD N NA OA"

# ---- copy the binaries and the receptor ------------------------------
for b in ("autodock4", "autogrid4"):
    shutil.copy2(os.path.join(SRC, b), os.path.join(AD4, b))
shutil.copy2(os.path.join(DOCK, "receptor.pdbqt"), os.path.join(AD4, "receptor.pdbqt"))

mini = pd.read_csv(os.path.join(DOCK, "docking_set_minimal.csv"))
lig_ids = list(mini["lig_id"])

for lid in lig_ids:
    shutil.copy2(os.path.join(DOCK, "ligands", lid + ".pdbqt"),
                 os.path.join(AD4, lid + ".pdbqt"))

# ---- grid parameter file --------------------------------------------
maps = " ".join(f"map receptor.{t}.map" for t in LIGAND_TYPES.split())
gpf = f"""npts {NPTS[0]} {NPTS[1]} {NPTS[2]}
gridfld receptor.maps.fld
spacing {SPACING}
receptor_types {RECEPTOR_TYPES}
ligand_types {LIGAND_TYPES}
receptor receptor.pdbqt
gridcenter {CENTER[0]} {CENTER[1]} {CENTER[2]}
smooth 0.5
"""
for t in LIGAND_TYPES.split():
    gpf += f"map receptor.{t}.map\n"
gpf += "elecmap receptor.e.map\ndsolvmap receptor.d.map\ndielectric -0.1465\n"
open(os.path.join(AD4, "grid.gpf"), "w", newline="\n").write(gpf)

# ---- one docking parameter file per ligand ---------------------------
DPF_TEMPLATE = """autodock_parameter_version 4.2
outlev 1
intelec
seed pid time
ligand_types {ltypes}
fld receptor.maps.fld
{maplines}
elecmap receptor.e.map
desolvmap receptor.d.map
move {lig}.pdbqt
about 0.0 0.0 0.0
tran0 random
quaternion0 random
dihe0 random
rmstol 2.0
extnrg 1000.0
e0max 0.0 10000
ga_pop_size 150
ga_num_evals 2500000
ga_num_generations 27000
ga_elitism 1
ga_mutation_rate 0.02
ga_crossover_rate 0.8
ga_window_size 10
ga_cauchy_alpha 0.0
ga_cauchy_beta 1.0
set_ga
sw_max_its 300
sw_max_succ 4
sw_max_fail 4
sw_rho 1.0
sw_lb_rho 0.01
ls_search_freq 0.06
set_psw1
unbound_model bound
ga_run 10
analysis
"""
maplines = "\n".join(f"map receptor.{t}.map" for t in LIGAND_TYPES.split())

for lid in lig_ids:
    # centre the ligand search on the grid centre
    dpf = DPF_TEMPLATE.format(ltypes=LIGAND_TYPES, maplines=maplines, lig=lid)
    dpf = dpf.replace("about 0.0 0.0 0.0",
                      f"about {CENTER[0]} {CENTER[1]} {CENTER[2]}")
    open(os.path.join(AD4, lid + ".dpf"), "w", newline="\n").write(dpf)

# ---- runner script ---------------------------------------------------
runner = """#!/bin/bash
set -e
cd "$(dirname "$0")"
chmod +x autodock4 autogrid4

echo "=== AUTOGRID: building affinity maps ==="
./autogrid4 -p grid.gpf -l grid.glg
echo "maps built"
echo

LIGANDS="%s"
n=0
total=$(echo $LIGANDS | wc -w)

for L in $LIGANDS; do
    n=$((n+1))
    if [ -s "dlg/$L.dlg" ]; then
        echo "[$n/$total] $L - already done, skipping"
        continue
    fi
    echo "[$n/$total] docking $L ..."
    ./autodock4 -p $L.dpf -l dlg/$L.dlg
done

echo
echo "=== ALL DOCKING COMPLETE ==="
""" % " ".join(lig_ids)

open(os.path.join(AD4, "run_ad4.sh"), "w", newline="\n").write(runner)

print("AutoDock4 setup written to:", AD4)
print("ligands        :", len(lig_ids))
print("grid centre    :", CENTER)
print("grid points    :", NPTS, "at", SPACING, "A  =>",
      round(NPTS[0] * SPACING, 1), "A per side")
print("receptor types :", RECEPTOR_TYPES)
print("ligand types   :", LIGAND_TYPES)
print("\nrun with:")
print('  wsl -d Ubuntu-22.04 -- bash /mnt/d/DES/CADD/Project/ad4/run_ad4.sh')
