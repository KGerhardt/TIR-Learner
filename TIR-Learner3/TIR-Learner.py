#!/usr/app/env python3
# -*- coding: utf-8 -*-
# Tianyu (Sky) Lu (tianyu@lu.fm)
# 2025-03-29

import os
import sys

sys.path.insert(0, f"{os.path.dirname(__file__)}/app")

# Use noqa: E402 to suppress "PEP8: E402" --------------------------------------------↴
import argparse                                                                 # noqa: E402
import shutil                                                                   # noqa: E402
from typing import List, Tuple, Optional                                        # noqa: E402

from app.main import TIRLearner                                                 # noqa: E402
from app.shared import DEFAULT_ALLOCATED_PROCESSORS, SKIP_TIRVISH, SKIP_GRF     # noqa: E402

VERSION = "v3.0.7"
INFO = "by Tianyu (Sky) Lu (tianyu@lu.fm) released under GPLv3"


def _process_additional_args(additional_args: Optional[List[str]]) -> Optional[Tuple[str, ...]]:
    if not additional_args:
        return None
    # processed_additional_args = tuple(map(str.upper, additional_args))
    processed_additional_args = tuple(map(lambda x: str(x.upper()), additional_args))
    if (SKIP_TIRVISH in processed_additional_args) and (SKIP_GRF in processed_additional_args):
        raise SystemExit("[ERROR] \"skip_tirvish\" and \"skip_grf\" cannot be specified at the same time!")
    return processed_additional_args


def main():
    """Main function to handle command line arguments and execute the program."""
    parser = argparse.ArgumentParser(prog="TIR-Learner",
                                     description="TIR-Learner is an ensemble pipeline for Terminal Inverted Repeat "
                                                 "(TIR) transposable elements annotation in eukaryotic genomes")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION} {INFO}")

    parser.add_argument("-f", "--genome_file", help="Genome file in fasta format",
                        type=str, required=True)
    parser.add_argument("-n", "--genome_name", help="Genome name (Optional)",
                        type=str, default="TIR-Learner")
    parser.add_argument("-s", "--species", help="One of the following: \"maize\", \"rice\" or \"others\"",
                        type=str, required=True)
    parser.add_argument("-l", "--length", help="Max length of TIR (Optional)", type=int, default=5000)
    parser.add_argument("-p", "--processors", "-t", "--cpu",
                        help="Number of processors allowed (Optional)", type=int, default=DEFAULT_ALLOCATED_PROCESSORS)
    # -t means --thread, however multithreading is abandoned, so it's only for downward compatibility
    # TODO add py and gnup two parallel execution mode, also add more detailed help info
    parser.add_argument("-m", "--mode", help=("Parallel execution mode, one of the following: \"py\" "
                                              "and \"gnup\" (Optional)"), type=str, default="py")
    parser.add_argument("-w", "--working_dir", help="The path to the working directory (Optional). "
                                                    "An isolated sandbox directory for storing all the temporary files "
                                                    "will be created in the working directory. This sandbox directory "
                                                    "will only persist during the program execution. DO NOT TOUCH "
                                                    "THE SANDBOX DIRECTORY IF IT IS NOT FOR DEBUGGING!",
                        type=str, default=None)
    parser.add_argument("-o", "--output_dir", help="Output directory (Optional)", type=str, default=None)
    parser.add_argument("-c", "--checkpoint_dir", help="The path to the checkpoint directory (Optional). "
                                                       "If not specified, the program will automatically search for it "
                                                       "in the genome file directory and the output directory.",
                        type=str, nargs='?', const="auto", default=None)
    parser.add_argument("--verbose", help="Verbose mode (Optional). "
                                          "Will show interactive progress bar and more execution details.",
                        action="store_true")
    parser.add_argument("-d", "--debug", help="Debug mode (Optional). If activated, data for all "
                                              "completed steps will be stored in the checkpoint file. Meanwhile, "
                                              "the temporary files in the working directory will also be kept.",
                        action="store_true")
    parser.add_argument("--grf_path", help="Path to GRF program (Optional)",
                        type=str, default=os.path.dirname(shutil.which("grf-main")))
    parser.add_argument("--gt_path", help="Path to genometools program (Optional)",
                        type=str, default=os.path.dirname(shutil.which("gt")))
    parser.add_argument("-a", "--additional_args", help="Additional arguments (Optional). "
                                                        "See documentation for more details.",
                        type=str, nargs="+")
    # see prog_const for what additional args are acceptable

    parsed_args = parser.parse_args()

    genome_file: str = parsed_args.genome_file
    output_dir: Optional[str] = parsed_args.output_dir
    if not output_dir:
        output_dir = os.path.dirname(genome_file)
    genome_file = os.path.abspath(genome_file)
    output_dir = os.path.abspath(output_dir)

    checkpoint_input: Optional[str] = parsed_args.checkpoint_dir
    if checkpoint_input and checkpoint_input != "auto":
        checkpoint_input = os.path.abspath(checkpoint_input)

    GRF_path: str = os.path.abspath(parsed_args.grf_path.replace('"', ""))
    gt_path: str = os.path.abspath(parsed_args.gt_path.replace('"', ""))

    additional_args: Optional[Tuple[str, ...]] = _process_additional_args(parsed_args.additional_args)
    if additional_args:
        print(f"[INFO] Additional args: {additional_args} captured.")

    TIRLearner_instance = TIRLearner(
        genome_file_path=genome_file,
        genome_name=parsed_args.genome_name,
        species=parsed_args.species,
        TIR_length=parsed_args.length,
        processors=parsed_args.processors,
        para_mode=parsed_args.mode,
        working_dir_path=parsed_args.working_dir,
        output_dir_path=output_dir,
        checkpoint_dir_input_path=checkpoint_input,
        flag_verbose=parsed_args.verbose,
        flag_debug=parsed_args.debug,
        GRF_path=GRF_path,
        gt_path=gt_path,
        additional_args=additional_args
    )
    TIRLearner_instance.execute()

if __name__ == "__main__":
    main()
