#!/bin/bash
set -e
cd "$(dirname "$0")"
chmod +x autodock4 autogrid4

echo "=== AUTOGRID: building affinity maps ==="
./autogrid4 -p grid.gpf -l grid.glg
echo "maps built"
echo

LIGANDS="E0001 E0002 E0003 E0004 E0005 E0006 E0007 E0008 E0009 E0010 E0011 E0012 E0013 E0014 L0000 L0001 L0002 L0003 L0004 L0005"
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
