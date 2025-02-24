import os
import warnings

os.environ["KERAS_BACKEND"] = "torch"  # use pytorch as keras backend
os.environ["KMP_WARNINGS"] = '0'  # mute all OpenMP warnings. #shujun
warnings.filterwarnings("ignore", category=UserWarning)  # mute keras warning

# Use if True to suppress the PEP8: E402 warning
if True:  # noqa: E402
    import gc
    import datetime
    import json
    import math
    import multiprocessing as mp
    import psutil
    import regex as re
    import shutil
    import subprocess
    import tempfile
    import time
    from typing import Union, Optional

    import numpy as np
    import pandas as pd
    import swifter

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    # Attention: sklearn does not automatically import its subpackages
    from sklearn.preprocessing import LabelEncoder

    import torch
    import keras

# Acceptable additional args
CHECKPOINT_OFF = "CHECKPOINT_OFF"
NO_PARALLEL = "NO_PARALLEL"
SKIP_TIRVISH = "SKIP_TIRVISH"
SKIP_GRF = "SKIP_GRF"

FILE_NAME_SPLITER = "-+-"
TERMINAL_SPLITER = '#'

TIR_SUPERFAMILIES = ("DTA", "DTC", "DTH", "DTM", "DTT")

CNN_MODEL_DIR_PATH = "./cnn0912/cnn0912.keras"
SANDBOX_DIR_NAME = "[DONT_ALTER]TIR-Learner_sandbox"
SPLIT_FASTA_TAG = "SplitFasta"

PROGRAM_ROOT_DIR_PATH = os.path.abspath(str(os.path.dirname(os.path.dirname(__file__))))
CNN_MODEL_DIR_PATH = os.path.join(PROGRAM_ROOT_DIR_PATH, CNN_MODEL_DIR_PATH)

REFLIB_DIR_NAME = "RefLib"
REFLIB_AVAILABLE_SPECIES = ("rice", "maize")
REFLIB_FILE_DICT = {species: [f"{species}_{TIR_type}_RefLib" for TIR_type in TIR_SUPERFAMILIES]
                    for species in REFLIB_AVAILABLE_SPECIES}
REFLIB_DIR_PATH = os.path.join(PROGRAM_ROOT_DIR_PATH, REFLIB_DIR_NAME)

DEFAULT_ALLOCATED_PROCESSORS = os.cpu_count() - 2 if os.cpu_count() > 2 else 1

TIRVISH_SPLIT_SEQ_LEN = 5 * (10 ** 6)  # 5 mb
TIRVISH_OVERLAP_SEQ_LEN = 50 * (10 ** 3)  # 50 kb

# TODO only for debug
# TIRvish_split_seq_len = 10 * (10 ** 3)
# TIRvish_overlap_seq_len = 10 * (10 ** 3)

SHORT_SEQ_LEN = 2000
# GENERAL_SPLIT_NUM_THRESHOLD = 5
# MIX_SPLIT_PERCENT_THRESHOLD = 0.05
# MIX_SHORT_SEQ_PROCESS_NUM = 2

DEFAULT_TERMINAL_WIDTH = 120
try:
    TERMINAL_WIDTH = os.get_terminal_size().columns
    if TERMINAL_WIDTH > DEFAULT_TERMINAL_WIDTH:
        TERMINAL_WIDTH = DEFAULT_TERMINAL_WIDTH
except OSError:
    TERMINAL_WIDTH = DEFAULT_TERMINAL_WIDTH


def terminal_print(string: str = "", spliter: str = TERMINAL_SPLITER, **kwargs):
    if string != "":
        num_spliter = (TERMINAL_WIDTH - (2 + len(string)) - 1) // 2
        print(f"{spliter * num_spliter} {string} {spliter * num_spliter}", **kwargs)
    else:
        print(spliter * (TERMINAL_WIDTH - 1), **kwargs)
