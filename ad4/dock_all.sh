#!/bin/bash
# Run AutoDock4 on every ligand. Maps must already exist (autogrid4 step).
# Skips ligands whose .dlg is already present, so it can be re-run safely.

cd "$(dirname "$0")"
mkdir -p dlg

LIGANDS=$(ls *.pdbqt | grep -v receptor | sed 's/\.pdbqt$//' | sort)
total=$(echo "$LIGANDS" | wc -l)
n=0

for L in $LIGANDS; do
    n=$((n+1))
    if [ -s "dlg/$L.dlg" ] && grep -q "CLUSTERING HISTOGRAM" "dlg/$L.dlg" 2>/dev/null; then
        echo "[$n/$total] $L - already done"
        continue
    fi
    echo "[$n/$total] docking $L ..."
    ./autodock4 -p "$L.dpf" -l "dlg/$L.dlg" > /dev/null 2>&1
    best=$(grep -E '^ +1 \|' "dlg/$L.dlg" | head -1 | awk -F'|' '{print $2}' | tr -d ' ')
    echo "        best cluster energy: $best kcal/mol"
done

echo
echo "=== AUTODOCK4 COMPLETE ==="
