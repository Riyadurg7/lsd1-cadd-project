"""
Local AutoDock Vina batch docking for the LSD1 project.

Usage:
    python run_docking.py                 # full run, exhaustiveness 8
    python run_docking.py --exh 4         # faster
    python run_docking.py --limit 20      # quick test on 20 ligands

Resumes automatically: ligands already scored in docking_results.csv are skipped.
"""
import argparse, csv, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

VINA     = os.path.join(HERE, "vina.exe")
RECEPTOR = "receptor.pdbqt"
LIGDIR   = "ligands"
OUTDIR   = "dockout"
RESULTS  = "docking_results.csv"
SETFILE  = "docking_set_prepared.csv"

ap = argparse.ArgumentParser()
ap.add_argument("--exh",   type=int, default=8,  help="exhaustiveness")
ap.add_argument("--cpu",   type=int, default=6,  help="CPU threads for Vina")
ap.add_argument("--chunk", type=int, default=25, help="ligands per invocation")
ap.add_argument("--limit", type=int, default=0,  help="stop after N ligands (0 = all)")
ap.add_argument("--box",   type=float, default=22.0, help="box edge in Angstrom")
ap.add_argument("--set",     dest="setfile", default="docking_set_prepared.csv",
                help="CSV listing the ligands to dock")
ap.add_argument("--results", default="docking_results.csv",
                help="CSV to append scores to")
args = ap.parse_args()

SETFILE = args.setfile
RESULTS = args.results


def box_from_crystal_ligand(path="crystal_ligand.pdb"):
    """Centre the search box on the crystallographic ligand."""
    xs, ys, zs = [], [], []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
    if not xs:
        sys.exit("ERROR: no atoms found in " + path)
    return [round(sum(v) / len(v), 3) for v in (xs, ys, zs)]


def best_score(path):
    try:
        for line in open(path):
            if line.startswith("REMARK VINA RESULT:"):
                return float(line.split()[3])
    except FileNotFoundError:
        pass
    return None


def load_rows():
    import csv as _csv
    with open(SETFILE, newline="", encoding="utf-8") as fh:
        return list(_csv.DictReader(fh))


def already_scored():
    if not os.path.exists(RESULTS):
        return set()
    out = set()
    with open(RESULTS, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            s = (row.get("score") or "").strip().lower()
            if s and s not in ("nan", "none", ""):
                out.add(row["lig_id"])
    return out


def main():
    for p in (VINA, RECEPTOR, SETFILE):
        if not os.path.exists(p):
            sys.exit("ERROR: missing " + p)

    os.makedirs(OUTDIR, exist_ok=True)
    centre = box_from_crystal_ligand()
    box = [args.box] * 3

    rows = load_rows()
    done = already_scored()
    todo = [r for r in rows
            if r["lig_id"] not in done
            and os.path.exists(os.path.join(LIGDIR, r["lig_id"] + ".pdbqt"))]
    if args.limit:
        todo = todo[:args.limit]

    print(f"receptor      : {RECEPTOR}")
    print(f"box centre    : {centre}")
    print(f"box size      : {box}")
    print(f"exhaustiveness: {args.exh}   cpu: {args.cpu}   chunk: {args.chunk}")
    print(f"already scored: {len(done)}")
    print(f"to dock       : {len(todo)}\n", flush=True)
    if not todo:
        print("nothing to do")
        return

    write_header = not os.path.exists(RESULTS)
    t0 = time.time()
    failures = 0

    with open(RESULTS, "a", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        if write_header:
            w.writerow(["lig_id", "group", "pIC50", "smiles", "score"])

        for start in range(0, len(todo), args.chunk):
            block = todo[start:start + args.chunk]
            paths = [os.path.join(LIGDIR, r["lig_id"] + ".pdbqt") for r in block]

            cmd = [VINA,
                   "--receptor", RECEPTOR,
                   "--batch", *paths,
                   "--dir", OUTDIR,
                   "--center_x", str(centre[0]),
                   "--center_y", str(centre[1]),
                   "--center_z", str(centre[2]),
                   "--size_x", str(box[0]),
                   "--size_y", str(box[1]),
                   "--size_z", str(box[2]),
                   "--exhaustiveness", str(args.exh),
                   "--num_modes", "5",
                   "--cpu", str(args.cpu),
                   "--seed", "42"]

            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                print("VINA FAILED (exit %d)" % proc.returncode)
                print(proc.stdout[-1500:])
                print(proc.stderr[-1500:])
                return

            for r in block:
                lid = r["lig_id"]
                score = best_score(os.path.join(OUTDIR, lid + "_out.pdbqt"))
                if score is None:
                    failures += 1
                w.writerow([lid, r["group"], r["pIC50"], r["clean_smiles"],
                            score if score is not None else "nan"])
            fh.flush()

            n = start + len(block)
            el = time.time() - t0
            eta = (len(todo) - n) * el / n / 60
            print(f"  {n}/{len(todo)}   {el/n:5.1f}s/ligand   ETA {eta:4.0f} min",
                  flush=True)

    print(f"\ndone - {failures} ligands produced no score")
    print("results:", os.path.join(HERE, RESULTS))


if __name__ == "__main__":
    main()
