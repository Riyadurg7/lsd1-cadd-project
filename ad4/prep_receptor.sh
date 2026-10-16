#!/bin/bash
# Prepare the LSD1 receptor for AutoDock4 using AutoDockTools.
#
# IMPORTANT: prepare_receptor4.py defaults to -U nphs_lps_waters_nonstdres,
# and "nonstdres" would delete FAD - the cofactor that sits at the bottom of
# the binding site. "nonstdres" is therefore omitted below, and FAD is
# verified in the output.

set -e

MGL=/mnt/d/DES/CADD/mgltools_x86_64Linux2_1.5.7
export MGL_ROOT="$MGL"
export PYTHONPATH="$MGL/MGLToolsPckgs"
export LD_LIBRARY_PATH="$MGL/lib:$LD_LIBRARY_PATH"

AD4=/mnt/d/DES/CADD/Project/ad4
cd "$AD4"

chmod +x "$MGL/bin/pythonsh" 2>/dev/null || true

PY="$MGL/bin/python2.7"
PREP="$MGL/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"

echo "=== preparing receptor (keeping FAD) ==="
"$PY" "$PREP" \
    -r receptor.pdb \
    -o receptor.pdbqt \
    -A checkhydrogens \
    -U nphs_lps_waters

echo
echo "=== verification ==="
echo -n "FAD atoms in output: "
grep -c "FAD" receptor.pdbqt || echo 0

echo -n "non-zero partial charges: "
grep -E "^(ATOM|HETATM)" receptor.pdbqt | awk '{print $(NF-1)}' | grep -vc "^+0.000$" || echo 0

echo -n "total atoms: "
grep -cE "^(ATOM|HETATM)" receptor.pdbqt

echo
echo "sample lines:"
grep -E "^(ATOM|HETATM)" receptor.pdbqt | head -3
echo "..."
grep "FAD" receptor.pdbqt | head -3
