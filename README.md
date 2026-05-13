# TIR-Learner v4

[![license](https://img.shields.io/github/license/KGerhardt/TIR-Learner.svg)](https://github.com/KGerhardt/TIR-Learner/blob/main/LICENSE)
[![bioconda platform](https://anaconda.org/bioconda/tir-learner/badges/platforms.svg)](https://anaconda.org/bioconda/tir-learner)
[![bioconda version](https://anaconda.org/bioconda/tir-learner/badges/version.svg)](https://anaconda.org/bioconda/tir-learner)
[![bioconda downloads](https://anaconda.org/bioconda/tir-learner/badges/downloads.svg)](https://anaconda.org/bioconda/tir-learner)

TIR-Learner is an ensemble pipeline for Terminal Inverted Repeat (TIR) transposable element annotation in eukaryotic genomes.

Version 4 is a complete, end-to-end rewrite of the program. It retains the logical workflow, the external software calls (BLAST, GRF, TIRvish, Keras), and the core CNN model from v3, but every supporting layer — genome handling, candidate processing, parallelization, intermediate I/O, and post-processing — has been reimplemented from scratch. 

Results in TIR-Learner v4 are essentially identical to those in v3, excepting a bugfix to correctly localize sequence retrievals from TIRvish, edge case inclusions/exclusions due to v4 correctly calculating TIR percent identity, and small changes introduced through (slight) indel allowance in TIR, TSD sequences. 

TIR-Learner v4 is about 3.3x faster than v3 and uses vastly less RAM, achieving <= 2GB RAM per thread usage. Processing genomes of any size is possible with v4.

## Table of Contents

- [Background](#background)
- [What Changed in Version 4](#what-changed-in-version-4)
- [Installation](#installation)
- [Usage](#usage)
- [Program Workflow](#program-workflow)
- [Output Files](#output-files)
- [Algorithm Details](#algorithm-details)
- [Citation](#citation)
- [Credits](#credits)
- [License](#license)

## Background

Transposable elements (TEs) are DNA sequences that can move within a genome. Terminal Inverted Repeat (TIR) transposons are a class of TE characterized by inverted repeat sequences at their ends. These TIRs act as bookends that mark the boundaries of the element and allow transposase enzymes to recognize and mobilize it.

TIR-Learner combines several approaches to identify and classify TIR transposons:
- Homology-based detection using curated reference libraries (rice and maize)
- De novo structural identification with TIRvish and GRF
- Convolutional neural network (CNN) classification into TIR superfamilies
- TSD/TIR pattern validation against superfamily-specific rules

## What Changed in Version 4

Version 4 keeps the v3 logical workflow (genome → TIRvish/GRF → optional homology → CNN → TSD/TIR validation → final annotation) and the trained CNN model (`cnn0912`) intact. Everything else has been rewritten:

### Functional alterations

- **Bugfixes.** V3 has a handful of sequence-retrieval and edge-case GRF/TIRvish parsing bugs that affect a minority of sequences. These are fixed in V4.
- **Improved TSD search behavior.** V3 used an inefficient TSD and TIR checking module that only accepted ungapped TIR/TSD sequences and prematurely terminated TIR similarity calculations after only 10 matching bases. V4 uses a more sensitive, indel-tolerant approach and correctly calculates TIR identity over the full span.
- **Transparent, persistent outputs.** All major computational outputs are kept in compact JSON summary files.

### Structural improvements

- **New genome splitter.** The input genome is chunked with overlap-aware splitting and indexed with `pyfastx`, replacing the previous splitter and the per-step re-reads of the genome.
- **Direct, in-memory candidate flow.** Candidates from TIRvish and GRF are emitted as compact JSON records and consumed directly by the CNN step, removing the v3 intermediate FASTA shuffling and most disk round-trips.
- **Native multiprocessing.** Parallelism is handled with the Python `multiprocessing` standard library, with TIRvish and GRF run per-chunk in parallel and CNN inference distributed across worker processes. The `swifter` and `keras-tuner` dependencies and the v3 `pyboost`/`pystrict`/`gnup` mode selector are gone.
- **Compressed checkpoints.** TIRvish and GRF outputs are written as gzipped JSON (`pigz` if available, `gzip` otherwise) and can be reloaded into a later run via `--existing_tirvish` / `--existing_grf` to skip recomputation.
- **Simplified CLI.** Most v3 flags have been removed in favor of a smaller, more direct argument set. Deprecated v3 flags (`-n`, `-m`, `-w`, `-c`, `--verbose`, `-d`, `--grf_path`, `--gt_path`, `-a`) are still accepted silently for backwards compatibility with the EDTA wrapper, but they have no effect (or, in the case of `-w`, alias `-o`).
- **CPU-first execution.** The CNN runs on CPU by default (`KERAS_BACKEND=torch`, `CUDA_VISIBLE_DEVICES` cleared), which avoids the GPU-build conda installation issues seen in v3.
- **Refreshed dependencies.** Keras 3 + PyTorch backend, NumPy 2.x, `pyfastx` for indexed genome access; the v3 `pandas`/`swifter` data path and `BioPy` have been removed.

### Performance comparison

TIR-Learner v3 suffered from a handful of issues that limited its usability on large genomes: 

* (3-1) The program took an aggressive, memory-forward method of filtering results using the Pandas package. While this approach was reasonably fast, memory copying in Python's parallel architecture and subtle design choices in TIR-Learner v3 often caused RAM consumption in the hundreds of GB for even reasonably sized genomes ~5 Gbp.
* (3-2) v3 used a very inefficient method of sequence retrieval to extract TIR candidates after processing. In smaller, more fragmented genomes (such as the Pacific shrimp genome used in TIR-Learner v3's testing), the scalability issues of the specific BioPy implementation were not readily apparent. Modern, complete chromosome genomes exposed this issue, and in most such cases, the outright majority of program runtime would be consumed in this computationally very simple process.

A typical eukaryotic genome with intact chromosomes would spend ~60% of its runtime in TIR-Learner v3 in these two steps, with TIRvish, GRF, and the CNN labelling collectively comprising the remaining 40%. While CPU utilization was good at ~70% @ 48 threads, the poor baseline performance of the underlying approaches made this figure misleading.

v4 solves these issues:

* (4-1) Fixed-size genome chunking separates original genome size from RAM use entirely.
* (4-2) More efficient filtering algorithms, a streaming filtering approach, and chunk-local filtering distribute work and ensure predictable RAM use/thread irrespective of genome size, TIR density.
* (4-3) Chunk-local sequence retrieval is ultra-fast, distributed, and requires little RAM.

The same eukaryotic genome processed with TIR-Learner v4 spends ~0.1% of its runtime in steps equivalent to (3-1) and (3-2). Although v4 is also more efficient in the GRF, TIRvish, and especially the CNN steps due to better load balancing and controlled resource allocation in the CNN stage, the baseline of work done by these tools is quite substantial. They account for nearly the entirety of the remaining ~30% of v4's runtime compared to v3. CPU utilization is also better, typically ~80% for shorter genomes and plateauing around ~90% for long ones.

## Installation

### Bioconda (recommended)

The bioconda `tir-learner` package now tracks v4. This is the preferred way to install — it pulls in BLAST, GenomeTools (TIRvish), GRF, PyTorch, and all other dependencies in one step.

*If you do not have conda/mamba installed, we recommend Miniforge: [conda-forge/miniforge](https://github.com/conda-forge/miniforge). Note that `-c conda-forge` MUST be specified before `-c bioconda` so that conda-forge takes priority (see [conda channel docs](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/channels.html#specifying-channels-when-installing-packages)).*

Install into a new environment:

```shell
mamba create -n TIRLearner_env -c conda-forge -c bioconda tir-learner
```

Or into an existing environment:

```shell
mamba install -c conda-forge -c bioconda tir-learner
```

After activation, the program is available on `PATH` as `TIR-Learner`.

### From source

Clone the repository and install dependencies via the bundled conda environment file:

```shell
git clone https://github.com/KGerhardt/TIR-Learner.git
cd TIR-Learner
mamba env create -f environment.yml -n TIRLearner_env
mamba activate TIRLearner_env
```
Run the program with:

```shell
python TIR-Learner4/TIR-Learner.py [options]
```

### Dependencies

- Python ≥ 3.10
- BLAST+
- GenomeTools (`gt`, provides TIRvish)
- GRF (Generic Repeat Finder)
- `pigz` (optional; falls back to `gzip`)
- Python packages:
  - `keras` ≥ 3
  - `pytorch`
  - `numpy`
  - `pyfastx`

See `environment.yml` for a full pinned dependency list.

## Usage

```shell
TIR-Learner.py -f GENOME_FILE [-s SPECIES] [-l LENGTH] [-p PROCESSORS] [-o OUTPUT_DIR]
               [--skip_tirvish] [--skip_grf]
               [--existing_tirvish JSON] [--existing_grf JSON]
```

### Required arguments

- `-f, -g, --genome_file`
  Genome file in FASTA format. Sequence names must be unique.

### Optional arguments

- `-s, --species`
  Reference library to use for the homology step:
  - `rice` — use the rice TIR reference library
  - `maize` — use the maize TIR reference library
  - `others` (or unset) — skip homology, run TIRvish/GRF + CNN only

- `-l, --length` (default: `5000`)
  Maximum length of TIR transposons to detect, in bp.

- `-p, -t, --cpu, --processors` (default: `1`)
  Number of parallel worker processes.

- `-o, --directory, --output_dir` (default: `./TIR_Learner_working_directory`)
  Working and output directory. Created if it does not exist. Final results are written to `<output_dir>/TIR-Learner-Result/`.

- `--skip_tirvish`
  Do not run TIRvish (and ignore any TIRvish results in post-processing).

- `--skip_grf`
  Do not run GRF (and ignore any GRF results in post-processing).

- `--existing_tirvish PATH`
  Path to a TIRvish JSON produced by a previous TIR-Learner run. The file is reused instead of re-running TIRvish.

- `--existing_grf PATH`
  Path to a GRF JSON produced by a previous TIR-Learner run. The file is reused instead of re-running GRF.

`--skip_tirvish` and `--skip_grf` cannot both be set; at least one de novo finder must run.

### Backwards compatibility

The following v3 flags are still accepted for compatibility with EDTA and existing scripts, but are silently ignored: `-n/--genome_name`, `-m/--mode`, `-c/--checkpoint_dir`, `--verbose`, `-d/--debug`, `--grf_path`, `--gt_path`, `-a/--additional_args`. `-w/--working_dir` is treated as an alias of `-o/--directory` when `-o` is at its default.

### Example

```bash
python TIR-Learner4/TIR-Learner.py -f ./test/test_4chr_1mb_rice.fa -s rice -p 8 -o ./tirlearner_out
```

## Program Workflow

<div style="text-align:center; width: 100%">
<br>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="./docs/TIR-Learner3_workflow_white.drawio.png">
  <img style="width: 100%; max-width: 5485px" alt="TIR-Learner workflow" src="./docs/TIR-Learner3_workflow_black.drawio.png">
</picture>
<br><br>
</div>

The high-level pipeline (unchanged from v3):

1. **Pre-scan and split.** The genome is indexed and split into overlapping chunks sized for downstream tools.
2. **De novo TIR detection.** TIRvish and GRF are run per chunk in parallel; their outputs are merged into compact JSON records.
3. **Optional homology (rice/maize only).** Candidates are checked against the species-specific reference library with BLAST. Hits are tagged; non-homologs are passed forward.
4. **CNN classification.** Each candidate is encoded and classified into one of the TIR superfamilies (`DTA`, `DTC`, `DTH`, `DTM`, `DTT`) or rejected as `NonTIR`.
5. **TSD/TIR validation.** Surviving candidates are checked against superfamily-specific TSD patterns and TIR conservation requirements.
6. **Post-processing.** Predictions are merged across chunks and emitted in two forms — a raw set and a deduplicated/overlap-resolved set.

## Output Files

All final outputs are written to `<output_dir>/TIR-Learner-Result/`.

### Annotations

- `TIR-Learner_FinalAnn.fa` — FASTA of all retained TIR predictions.
- `TIR-Learner_FinalAnn.gff3` — GFF3 with location and superfamily for each prediction.
- `TIR-Learner_FinalAnn_filter.fa` — FASTA after overlap resolution.
- `TIR-Learner_FinalAnn_filter.gff3` — GFF3 after overlap resolution.

### Reusable intermediates

These gzipped JSON files are checkpoints from the de novo finders and the post-CNN stage. They can be supplied to a later run via `--existing_tirvish` / `--existing_grf` to skip recomputation.

- `TIRVish_json.txt.gz` — raw TIRvish candidates.
- `GRF_json.txt.gz` — raw GRF candidates.
- `TIRVish_json_no_homologs.txt.gz`, `GRF_json_no_homologs.txt.gz` — candidates remaining after the homology step (rice/maize only).
- `post_CNN_TIRVish_json.txt.gz`, `post_CNN_GRF_json.txt.gz` — candidates that passed CNN classification and TSD/TIR validation.
- `Module1_homology_hits_against_genome_*` — BLAST homology hits, when the homology step ran.

## Algorithm Details

### TIR detection ranges

- Minimum TIR length: 10 bp
- Maximum TIR length: 1000 bp
- Maximum TIR-to-TIR distance: 5000 bp (configurable via `-l`)
- Minimum TIR similarity: 80%

### TSD validation

Superfamily-specific TSD patterns:

| Superfamily | TSD |
|-------------|-----|
| DTA | 8 bp |
| DTC | 2–3 bp |
| DTH | 3 bp (TWA) |
| DTM | 7–10 bp |
| DTT | 2 bp (TA) |

### CNN classification

- Trained model: `cnn0912` (Keras 3 / PyTorch backend), shipped under `TIR-Learner4/cnn0912/`.
- Input: encoded sequence fragments drawn from the TIR ends and adjacent regions of each candidate.
- Output: one of `DTA`, `DTC`, `DTH`, `DTM`, `DTT`, or `NonTIR`.

## Citation

TIR-Learner v4 has not yet been published and will likely appear bundled alongside the release of a new EDTA version.

The manuscript describing TIR-Learner v3 is in preparation. Slides:

[TIR-Learner v3: New generation TE annotation program for identifying TIRs](https://doi.org/10.6084/m9.figshare.26082574.v1)

Previous publications:

- TIR-Learner v2 (as part of EDTA v1):
  [Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline](https://doi.org/10.1186/s13059-019-1905-y)
- TIR-Learner v1:
  [TIR-Learner, a New Ensemble Method for TIR Transposable Element Annotation, Provides Evidence for Abundant New Transposable Elements in the Maize Genome](https://doi.org/10.1016/j.molp.2019.02.008)

## Credits

### Authors

- **Version 4 (rewrite):** [Kenji Gerhardt](https://github.com/KGerhardt)
- **Version 3:** [Tianyu (Sky) Lu](https://github.com/lutianyu2001)
- **Versions 1 and 2:** [Weijia Su](https://github.com/WeijiaSu) and [Shujun Ou](https://github.com/oushujun)

The v4 rewrite preserves the v3 workflow design and the trained CNN model; the underlying implementation is otherwise new.

### Acknowledgments

Development has been supported by:
- [The Ou Lab at The Ohio State University](https://www.ou-lab.org/)
- [The Ohio Supercomputer Center](https://www.osc.edu/)
- [OSU Undergraduate Research Access Innovation Seed Grant](https://ugresearch.osu.edu/faculty/funding-and-grants)

## License

This project is licensed under the GPL-3.0 License — see the [LICENSE](LICENSE) file for details.
