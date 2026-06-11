# TIRvish port — gold-set oracle fixture

Fixed test fixture for the faithful Rust port of `gt tirvish` (see the project
memory `tirvish-rust-port-plan`). The contract: the Rust tool must reproduce
`gt tirvish`'s **raw prediction multiset** on identical input. These chunks +
gold sets are that input and that expected output, committed so the acceptance
test is reproducible without re-running the ~9 min/chunk `gt tirvish` baseline.

## Contents (committed)

- `chunks/chunk{0..3}.fa` — fixed input: each a multi-FASTA of consecutive real
  Pacific white shrimp contigs (~5 Mb, records named `<contig>;;0`, uppercased),
  mirroring how the TIR-Learner splitter batches small contigs into a chunk.
- `chunks/manifest.tsv` — per-chunk contig count and bp.
- `gold/chunk{0..3}.tirvish.gff` — raw `gt tirvish` GFF3 output (the ground truth).
- `gold/chunk{0..3}.gold.tsv` — parsed gold set, one row per prediction:
  `seqid start stop tir1 tir2 tsd1 tsd2 sim`. First 7 cols = comparison KEY
  (full element coords incl. TSDs, 1-based; TIR/TSD lengths). `sim` = gt's
  tir_similarity, a diagnostic (not part of the key).

Not committed (see `.gitignore`): `idx/` (regenerable suffixerator indices),
`*.err`.

## Pipeline

```
make_chunks.py      # build chunks/chunk*.fa from ../pacific_white_shrimp.full.fa
run_gt_tirvish.sh   # gt suffixerator -mirrored + gt tirvish (EXACT pipeline flags) -> gold/*.tirvish.gff
parse_tirvish.py    # gold/*.tirvish.gff -> gold/*.gold.tsv (raw tirvish predictions, NO TA/N%/WFA filtering)
compare.py          # acceptance test: diff a candidate TSV against a gold TSV
```

`gt` flags are pinned to exactly what `tirvish_new.one_tirvish` runs:
`suffixerator -tis -suf -lcp -des -ssp -sds -dna -mirrored`, then
`tirvish -seed 20 -mintirlen 10 -maxtirlen 1000 -mintirdist 10 -maxtirdist 5000
-similar 80 -mintsd 2 -maxtsd 11 -vic 13 -seqids yes`.

## Using the oracle (once the Rust tool exists)

Emit the same TSV from the Rust tool on `chunks/chunkN.fa`, then:

```
python3 compare.py gold/chunkN.gold.tsv cand/chunkN.tsv --exact     # byte-exact coords+lengths
python3 compare.py gold/chunkN.gold.tsv cand/chunkN.tsv --tol 5     # boundary agreement within +/-5 bp
```

Reports recall / precision / missed (gold-only) / extra (cand-only).

## Baseline (gt tirvish wall time, the thing we're replacing)

| chunk | bp | contigs | gt tirvish wall | predictions |
|-------|-----|---------|-----------------|-------------|
| chunk0 | 5,248,098 | 16 | 573s (9m33s) | 270 |
| chunk1 | 5,431,767 | 10 | 695s | 267 |
| chunk2 | 5,056,285 | 16 | 500s | 295 |
| chunk3 | 5,702,366 | 14 | 482s | 321 |

Total: 1153 predictions over ~21 Mb. ~8-12 min / 5 Mb on repeat-rich shrimp
(chunk1, 10 large contigs, was the worst) = why full-genome tirvish is ~1.5 hr.
suffixerator is ~6 s of that; tirvish search is ~99% of the wall.
