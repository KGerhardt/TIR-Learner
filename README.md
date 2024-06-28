***---*** ⚠ ***Readme under construction ---***

# TIR-Learner

[![current release](https://img.shields.io/github/release/lutianyu2001/TIR-Learner.svg)](https://github.com/lutianyu2001/TIR-Learner/releases)
[![license](https://img.shields.io/github/license/lutianyu2001/TIR-Learner.svg)](https://github.com/lutianyu2001/TIR-Learner/blob/master/LICENSE)
[![run tests](https://github.com/lutianyu2001/TIR-Learner/actions/workflows/test.yml/badge.svg)](https://github.com/lutianyu2001/TIR-Learner/actions/workflows/test.yml)
[![bioconda platform](https://anaconda.org/bioconda/tir-learner/badges/platforms.svg)](https://anaconda.org/bioconda/tir-learner) 
[![bioconda version](https://anaconda.org/bioconda/tir-learner/badges/version.svg)](https://anaconda.org/bioconda/tir-learner) 
[![bioconda downloads](https://anaconda.org/bioconda/tir-learner/badges/downloads.svg)](https://anaconda.org/bioconda/tir-learner)

## Installation with Conda/Mamba

```bash
conda create -n TIR-Learner
conda activate TIR-Learner
mamba install -c conda-forge -c bioconda tir-learner
```

## Quick Test

```bash
TIR-Learner -f test/genome.fa -s others -t 2 -l 5000 -c ./ -o results
```
