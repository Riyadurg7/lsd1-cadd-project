"""
Enamine similarity search for the LSD1 project.

Same approach as the reference workflow (MACCS keys + Tanimoto, chunked over
the Enamine catalogue), with three differences:

  1. Salts are stripped - the Enamine file is full of entries like "Cl.N#CCN"
  2. Tiny fragments are skipped - they cannot resemble a drug-sized inhibitor
  3. Hits are collected from 0.70 upwards, not only above 0.90

Point 3 matters: a strict 0.90 cutoff can return nothing at all, which leaves
you with no ligands and no information. Collecting from 0.70 lets you report
counts at every threshold and still choose 0.90 if it yields enough.

Usage:
    python enamine_similarity.py                    # full run
    python enamine_similarity.py --limit 200000     # quick test
    python enamine_similarity.py --fp morgan        # Morgan instead of MACCS
"""
import argparse, csv, os, sys, time
from multiprocessing import Pool

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # project root

REF_FILE = os.path.join(HERE, "enamine", "LSD1_reference_compounds.csv")
ENA_FILE = os.path.join(HERE, "enamine", "Enamine_SMILES.csv")
OUT_FILE = os.path.join(HERE, "enamine", "Enamine_similarity_hits.csv")

# globals populated in each worker process
REF_FPS, REF_SMI, MIN_HEAVY, COLLECT, FPKIND, MORGAN = None, None, 12, 0.70, "maccs", None


def build_fp(mol):
    if FPKIND == "maccs":
        return MACCSkeys.GenMACCSKeys(mol)
    return MORGAN.GetFingerprint(mol)


def init_worker(ref_fps, ref_smi, min_heavy, collect, fpkind):
    global REF_FPS, REF_SMI, MIN_HEAVY, COLLECT, FPKIND, MORGAN
    REF_FPS, REF_SMI = ref_fps, ref_smi
    MIN_HEAVY, COLLECT, FPKIND = min_heavy, collect, fpkind
    MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def process_block(smiles_list):
    """Return [(enamine_smiles, best_similarity, reference_smiles), ...]"""
    hits = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
            if len(frags) > 1:
                mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
            if mol.GetNumHeavyAtoms() < MIN_HEAVY:
                continue
            fp = build_fp(mol)
            sims = DataStructs.BulkTanimotoSimilarity(fp, REF_FPS)
            best = max(sims)
            if best >= COLLECT:
                hits.append((smi, round(best, 4), REF_SMI[sims.index(best)]))
        except Exception:
            continue
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fp", choices=["maccs", "morgan"], default="maccs")
    ap.add_argument("--collect", type=float, default=0.70,
                    help="record every hit at or above this similarity")
    ap.add_argument("--min-heavy", type=int, default=12,
                    help="skip molecules with fewer heavy atoms than this")
    ap.add_argument("--refs", type=int, default=25,
                    help="how many reference actives to use")
    ap.add_argument("--cpu", type=int, default=6)
    ap.add_argument("--chunk", type=int, default=200000)
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after N catalogue compounds (0 = all)")
    args = ap.parse_args()

    for p in (REF_FILE, ENA_FILE):
        if not os.path.exists(p):
            sys.exit("ERROR: missing " + p)

    # ---- reference compounds -------------------------------------------
    ref = pd.read_csv(REF_FILE).head(args.refs)
    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ref_fps, ref_smi = [], []
    for smi in ref["SMILES"].dropna():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        ref_fps.append(MACCSkeys.GenMACCSKeys(m) if args.fp == "maccs"
                       else morgan.GetFingerprint(m))
        ref_smi.append(smi)

    print(f"fingerprint        : {args.fp.upper()}")
    print(f"reference actives  : {len(ref_fps)}")
    print(f"collecting hits >= : {args.collect}")
    print(f"min heavy atoms    : {args.min_heavy}")
    print(f"workers            : {args.cpu}\n", flush=True)

    if os.path.exists(OUT_FILE):
        os.remove(OUT_FILE)

    total, kept = 0, 0
    best_overall = 0.0
    t0 = time.time()

    pool = Pool(args.cpu, initializer=init_worker,
                initargs=(ref_fps, ref_smi, args.min_heavy,
                          args.collect, args.fp))

    with open(OUT_FILE, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Enamine_SMILES", "Tanimoto", "Reference_SMILES"])

        reader = pd.read_csv(ENA_FILE, chunksize=args.chunk, usecols=["SMILES"])
        for nchunk, df in enumerate(reader, start=1):
            smiles = df["SMILES"].dropna().tolist()

            per = max(1, len(smiles) // args.cpu + 1)
            blocks = [smiles[i:i + per] for i in range(0, len(smiles), per)]

            for res in pool.imap_unordered(process_block, blocks):
                for row in res:
                    w.writerow(row)
                    kept += 1
                    best_overall = max(best_overall, row[1])
            fh.flush()

            total += len(smiles)
            el = time.time() - t0
            print(f"chunk {nchunk:>3} | processed {total:>9,} | "
                  f"hits >= {args.collect}: {kept:>6,} | "
                  f"best so far {best_overall:.3f} | {el/60:5.1f} min",
                  flush=True)

            if args.limit and total >= args.limit:
                print("\nreached --limit, stopping early")
                break

    pool.close()
    pool.join()

    print("\n" + "=" * 58)
    print("SIMILARITY SEARCH COMPLETE")
    print("=" * 58)
    print(f"compounds processed : {total:,}")
    print(f"hits recorded       : {kept:,}")
    print(f"highest similarity  : {best_overall:.3f}")

    if kept:
        h = pd.read_csv(OUT_FILE)
        print("\nHITS BY THRESHOLD")
        for t in (0.95, 0.90, 0.85, 0.80, 0.75, 0.70):
            print(f"   >= {t:.2f} : {int((h['Tanimoto'] >= t).sum()):>7,}")
        print("\nTOP 10 HITS")
        print(h.nlargest(10, "Tanimoto")[["Tanimoto", "Enamine_SMILES"]]
               .to_string(index=False))
    print("\nsaved to:", OUT_FILE)


if __name__ == "__main__":
    main()
