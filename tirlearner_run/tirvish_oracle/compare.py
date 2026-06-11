#!/usr/bin/env python3
"""
Compare a candidate set against the gt tirvish gold set, keyed on the prediction
tuple. This is the acceptance test for the Rust TIRvish port: run it on the same
oracle chunks, emit the same TSV (seqid start stop tir1 tir2 tsd1 tsd2), and this
reports how the multisets agree.

Two match modes:
  --exact  : key = (seqid, start, stop, tir1, tir2, tsd1, tsd2)  [byte-for-byte]
  --tol N  : key = (seqid, round(start,N), round(stop,N)) only; boundaries within
             +/-N bp count as a match (default N=0 -> exact coords, lengths ignored)

Reports, per file pair and pooled: |gold|, |cand|, matched, gold-only (missed
recall), cand-only (false positives), and recall/precision. Exact coordinate
equality over millions of bases is not expected without full algorithmic
convergence; use --tol to see how close the boundaries land.

Usage:
  python3 compare.py gold.tsv cand.tsv [--exact | --tol N]
"""
import sys
from collections import Counter


def load(path):
    rows = []
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 7:
                continue
            seqid, start, stop, tir1, tir2, tsd1, tsd2 = p[:7]
            rows.append((seqid, int(start), int(stop), int(tir1), int(tir2), int(tsd1), int(tsd2)))
    return rows


def key_exact(r):
    return r  # full tuple


def key_tol(r, n):
    # bucket coords to +/-n; lengths ignored. (Bucketing, not interval-overlap:
    # cheap and good enough to gauge boundary agreement.)
    seqid, start, stop = r[0], r[1], r[2]
    if n <= 0:
        return (seqid, start, stop)
    return (seqid, round(start / (2 * n + 1)), round(stop / (2 * n + 1)))


def main():
    gold_path, cand_path = sys.argv[1], sys.argv[2]
    mode = "exact"
    tol = 0
    if "--tol" in sys.argv:
        mode = "tol"; tol = int(sys.argv[sys.argv.index("--tol") + 1])
    elif "--exact" in sys.argv:
        mode = "exact"

    gold = load(gold_path)
    cand = load(cand_path)
    kf = (lambda r: key_exact(r)) if mode == "exact" else (lambda r: key_tol(r, tol))

    gc, cc = Counter(map(kf, gold)), Counter(map(kf, cand))
    matched = sum((gc & cc).values())
    gold_only = sum((gc - cc).values())
    cand_only = sum((cc - gc).values())
    recall = matched / len(gold) if gold else float('nan')
    prec = matched / len(cand) if cand else float('nan')

    print(f"mode={mode}" + (f"(+/-{tol}bp)" if mode == "tol" else ""))
    print(f"  gold={len(gold)}  cand={len(cand)}  matched={matched}")
    print(f"  gold-only (missed) ={gold_only}   cand-only (extra) ={cand_only}")
    print(f"  recall={recall:.4f}  precision={prec:.4f}")
    return matched, gold_only, cand_only


if __name__ == "__main__":
    main()
