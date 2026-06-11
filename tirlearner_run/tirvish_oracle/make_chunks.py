#!/usr/bin/env python3
"""
Build a handful of fixed oracle chunks from the Pacific white shrimp genome.

Mirrors the TIR-Learner splitter closely enough for a faithful gt-tirvish test:
each chunk is a MULTI-FASTA of consecutive contigs (the real pipeline batches
small contigs into a ~chunk_size file), every record named `<contig>;;0` (the
`name;;offset` convention; offset 0 since we don't sub-split), sequence
uppercased (v4.07 genomeSplitter per-slice upper()). No concatenation across
contigs -> no spurious cross-boundary TIRs; gt tirvish processes each record
independently, exactly as in the pipeline.

The chunks are a *fixed test fixture*: gt tirvish defines the gold set on them,
and the future Rust tool runs on the identical files. So precise fidelity to the
splitter's packing is irrelevant -- only that it's real shrimp sequence.

Usage: python3 make_chunks.py [genome.fa] [n_chunks] [chunk_bp]
"""
import sys, os
import pyfastx

HERE = os.path.dirname(os.path.abspath(__file__))
GENOME = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "pacific_white_shrimp.full.fa")
N_CHUNKS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
CHUNK_BP = int(sys.argv[3]) if len(sys.argv) > 3 else 5_000_000

fa = pyfastx.Fasta(GENOME)
outdir = os.path.join(HERE, "chunks")

# Pack consecutive contigs (file order) into ~CHUNK_BP multi-FASTA chunks.
chunk_idx = 0
cur, cur_bp = [], 0
manifest = []
for i in range(len(fa)):
    rec = fa[i]
    cur.append((rec.name, rec.seq.upper()))
    cur_bp += len(rec)
    if cur_bp >= CHUNK_BP:
        path = os.path.join(outdir, f"chunk{chunk_idx}.fa")
        with open(path, "w") as fh:
            for name, seq in cur:
                fh.write(f">{name};;0\n")
                for j in range(0, len(seq), 80):
                    fh.write(seq[j:j+80] + "\n")
        manifest.append((chunk_idx, len(cur), cur_bp))
        print(f"chunk{chunk_idx}: {len(cur)} contigs, {cur_bp:,} bp -> {path}")
        chunk_idx += 1
        cur, cur_bp = [], 0
        if chunk_idx >= N_CHUNKS:
            break

with open(os.path.join(HERE, "chunks", "manifest.tsv"), "w") as fh:
    fh.write("chunk\tn_contigs\tbp\n")
    for c, n, bp in manifest:
        fh.write(f"chunk{c}\t{n}\t{bp}\n")
print(f"Wrote {len(manifest)} chunks.")
