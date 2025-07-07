#!/usr/app/env python3
# -*- coding: utf-8 -*-

import os
import warnings

os.environ["KERAS_BACKEND"] = "torch"  # use pytorch as keras backend
os.environ["KMP_WARNINGS"] = '0'  # mute all OpenMP warnings
warnings.filterwarnings("ignore", category=UserWarning)  # mute keras warning

# Use noqa to suppress false positive "Unused import statement" and "PEP8: E402" --↴
import gc                                                                       # noqa
import datetime                                                                 # noqa
import json                                                                     # noqa
import math                                                                     # noqa
import multiprocessing as mp                                                    # noqa
import psutil                                                                   # noqa
import regex as re                                                              # noqa
import shutil                                                                   # noqa
import subprocess                                                               # noqa
import tempfile                                                                 # noqa
import time                                                                     # noqa
from typing import Set, Tuple, List, Dict, Iterable, Sequence, Union, Optional  # noqa

import numpy as np                                                              # noqa
import pandas as pd                                                             # noqa
import swifter                                                                  # noqa

from Bio import SeqIO                                                           # noqa
from Bio.Seq import Seq                                                         # noqa
from Bio.SeqRecord import SeqRecord                                             # noqa

# Attention: sklearn does not automatically import its subpackages
from sklearn.preprocessing import LabelEncoder                                  # noqa

import torch                                                                    # noqa
import keras                                                                    # noqa

from collections import Counter

# ======================================================================================================================
# Constants

PROGRAM_ROOT_DIR_PATH: str = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

# Acceptable additional args
CHECKPOINT_OFF: str = "CHECKPOINT_OFF"
NO_PARALLEL: str = "NO_PARALLEL"
SKIP_TIRVISH: str = "SKIP_TIRVISH"
SKIP_GRF: str = "SKIP_GRF"

SPLITER: str = "-+-"
TERMINAL_SPLITER: str = '#'

ACCEPTED_NUCLEOTIDES: Set[str] = {'A', 'T', 'G', 'C', 'N'}
TIR_SUPERFAMILIES: Tuple[str, ...] = ("DTA", "DTC", "DTH", "DTM", "DTT")

SANDBOX_DIR_NAME: str = "[DONT_ALTER]TIR-Learner_sandbox"
SPLIT_FASTA_TAG: str = "SplitFasta"
RESULT_OUTPUT_DIR_NAME: str = "TIR-Learner-Result"
CHECKPOINT_DIR_NAME_PREFIX: str = "TIR-Learner_v3_checkpoint_"
PROCESSED_DE_NOVO_RESULT_FILE_NAME_FORMAT_STR: str = "{0}" + SPLITER + "processed_de_novo_result.fa"
PRINT_DECIMAL_SIGNIFICANT_FIGURES = 4

CNN_MODEL_DIR_REL_PATH: str = "./cnn0912/cnn0912.keras"
CNN_MODEL_DIR_ABS_PATH: str = os.path.join(PROGRAM_ROOT_DIR_PATH, CNN_MODEL_DIR_REL_PATH)

REFLIB_DIR_NAME: str = "RefLib"
REFLIB_AVAILABLE_SPECIES: Tuple[str, ...] = ("rice", "maize")
REFLIB_FILE_DICT: Dict[str, List[str]] = {species: [f"{species}_{TIR_type}_RefLib" for TIR_type in TIR_SUPERFAMILIES]
                                          for species in REFLIB_AVAILABLE_SPECIES}
REFLIB_DIR_PATH: str = os.path.join(PROGRAM_ROOT_DIR_PATH, REFLIB_DIR_NAME)

DEFAULT_ALLOCATED_PROCESSORS: int = os.cpu_count() - 2 if os.cpu_count() > 2 else 1

SHORT_SEQ_LEN: int = 2000

MP_SPLIT_SEQ_LEN: int = 5 * int(1e6)  # 5 Mbp
MP_OVERLAP_SEQ_LEN: int = 50 * int(1e3)  # 50 kbp
MP_SPLIT_SEQ_TAG: str = "_TIR-Learner_MP_split_"
MP_SPLIT_SEQ_ID_FORMAT_STR: str = "{0}" + MP_SPLIT_SEQ_TAG + "{1}of{2}"

# Debug only
# MP_SPLIT_SEQ_LEN = 10 * int(1e3)  # 10 kbp
# MP_OVERLAP_SEQ_LEN = 10 * int(1e3)  # 10 kbp

TE_TOLERANCE = 200

DEFAULT_TERMINAL_WIDTH: int = 120
try:
    TERMINAL_WIDTH = os.get_terminal_size().columns
    if TERMINAL_WIDTH > DEFAULT_TERMINAL_WIDTH:
        TERMINAL_WIDTH = DEFAULT_TERMINAL_WIDTH
except OSError:
    TERMINAL_WIDTH = DEFAULT_TERMINAL_WIDTH

# ======================================================================================================================
# Functions

def terminal_print(string: str = "", spliter: str = TERMINAL_SPLITER, **kwargs):
    if string != "":
        num_spliter = (TERMINAL_WIDTH - (2 + len(string)) - 1) // 2
        print(f"{spliter * num_spliter} {string} {spliter * num_spliter}", **kwargs)
    else:
        print(spliter * (TERMINAL_WIDTH - 1), **kwargs)
