[![Run tests](https://github.com/GallVp/TIR-Learner/actions/workflows/test.yml/badge.svg)](https://github.com/GallVp/TIR-Learner/actions/workflows/test.yml)

# TIR-Learner

## Run with Conda

```bash
conda env create --file environment.yml
conda activate tir-learner

./TIR-Learner3.0/TIR-Learner3.0.py \
    -f test/genome.fa \
    -s others \
    -t 2 \
    -l 5000 \
    -c ./ \
    -o results
```
