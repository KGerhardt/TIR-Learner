# import math
# import os
# import shutil
# import subprocess
# import pandas as pd
# import swifter  # ATTENTION: DO NOT REMOVE "swifter" EVEN IF IDE SHOWS IT IS NOT USED!
# import multiprocessing as mp
# from Bio import SeqIO
#
# import prog_const
# import os

from const import *

from run_TIRvish import process_fasta, retrieve_split_sequence_offset


def GRF(genome_file, genome_name, cpu_cores, TIR_length, GRF_path):
    GRF_bin_path = os.path.join(GRF_path, "grf-main")
    GRF_result_dir_name = f"{genome_name}_GRFmite"
    grf = (f"\"{GRF_bin_path}\" -i \"{genome_file}\" -o \"{GRF_result_dir_name}\" -c 1 -t {int(cpu_cores)} -p 20 "
           f"--min_space 10 --max_space {int(TIR_length)} --max_indel 0 --min_tr 10 "
           f"--min_spacer_len 10 --max_spacer_len {int(TIR_length)}")
    shell_filter = r" | grep -vE 'start:|end:|print:|calculate|^$'"

    # TODO debug only
    # print(grf)
    # print(shell_filter)
    subprocess.Popen(grf + shell_filter, shell=True).wait()


def GRF_mp(genome_file_path, genome_name, cpu_cores, TIR_length, GRF_path):
    GRF_working_dir = os.path.dirname(genome_file_path)
    genome_file_name = os.path.basename(genome_file_path)
    os.chdir(GRF_working_dir)
    GRF(genome_file_name, genome_name, cpu_cores, TIR_length, GRF_path)
    os.chdir("../")


def cpu_cores_allocation_GRF_boost(cpu_cores, num_seq):
    num_threads = int(cpu_cores * thread_core_ratio)
    num_process = int(math.sqrt(num_threads/16) + 1) * int(1 + num_seq/2)
    # num_process = floor((num_threads/16)^(1/2) + 1) * floor(1 + num_seq/2)
    return num_process, num_threads


def run_GRF_native(genome_file, genome_name, cpu_cores, TIR_length, flag_debug, GRF_path):
    print("  Step 1/2: Executing GRF in native mode\n")
    GRF(genome_file, genome_name, cpu_cores, TIR_length, GRF_path)
    print("  Step 2/2: Getting GRF result")
    return get_GRF_result_df_native(genome_name, flag_debug)


def run_GRF_boost(genome_file, genome_name, cpu_cores, TIR_length, flag_debug, GRF_path, fasta_files_path_list):
    os.makedirs(f"{splited_fasta_tag}_mp", exist_ok=True)  # TODO revise to include mode detection function
    os.chdir(f"./{splited_fasta_tag}_mp")

    print("  Step 1/3: Checking processed FASTA files")
    if (len(fasta_files_path_list) == 0 or
            any(not os.path.exists(f) or os.path.getsize(f) == 0 for f in fasta_files_path_list)):
        print("Processed FASTA files not found / invalid. Re-process FASTA files.")
        fasta_files_path_list.extend(process_fasta(genome_file, TIRvish_split_seq_len, TIRvish_overlap_seq_len))

    print("  Step 2/3: Executing GRF in boost mode\n")
    mp_args_list = [(file_path, genome_name, 2, TIR_length, GRF_path) for file_path in fasta_files_path_list]
    # TODO cpu_cores algorithm needed here
    with mp.Pool(cpu_cores) as pool:
        pool.starmap(GRF_mp, mp_args_list)
    print()

    print("  Step 3/3: Getting GRF result")
    return get_GRF_result_df_boost(fasta_files_path_list, genome_name, flag_debug,
                                   TIRvish_split_seq_len, TIRvish_overlap_seq_len)


# TODO serial to parallel
# def get_GRF_result_df_boost(fasta_files_path_list, genome_name, flag_debug, split_seq_len, overlap_seq_len):
#     GRF_result_dir_name = f"{genome_name}_GRFmite"
#     GRF_result_dir_list = [os.path.join(os.path.dirname(file), GRF_result_dir_name) for file in
#                            fasta_files_path_list]
#
#     df_list = []
#     for f in GRF_result_dir_list:
#         try:
#             df_data_dict = [{"id": rec.id, "seq": str(rec.seq), "len": len(rec)}
#                             for rec in SeqIO.parse(os.path.join(f, "candidate.fasta"), "fasta")]
#             df_in = pd.DataFrame(df_data_dict, columns=["id", "seq", "len"]).astype({"len": int})
#
#             id_pattern = r"^(\w+)_split_([\w.]+of\d+):(\d+):(\d+):(\w+):(\w+)$"
#
#             def revise_id(match):
#                 seqid, segment_position, start, end, TIR_pattern, TSD = match.groups()
#                 offset = retrieve_split_sequence_offset(segment_position, split_seq_len, overlap_seq_len)
#                 return f"{seqid}:{int(start) + offset}:{int(end) + offset}:{TIR_pattern}:{TSD}"
#
#             df_in["id"] = df_in["id"].str.replace(id_pattern, revise_id, regex=True)
#
#             df_list.append(df_in.copy())
#         except FileNotFoundError:
#             continue
#         if not flag_debug:
#             subprocess.Popen(["rm", "-rf", f])
#     return pd.concat(df_list).drop_duplicates(ignore_index=True).sort_values("id", ignore_index=True)


def get_GRF_result_df_boost(fasta_files_path_list, genome_name, flag_debug, split_seq_len, overlap_seq_len):
    GRF_result_dir_name = f"{genome_name}_GRFmite"
    GRF_result_dir_list = [os.path.join(os.path.dirname(file), GRF_result_dir_name) for file in
                           fasta_files_path_list]

    df_list = []
    for f in GRF_result_dir_list:
        try:
            df_data_dict = [{"id": rec.id, "seq": str(rec.seq), "len": len(rec)}
                            for rec in SeqIO.parse(os.path.join(f, "candidate.fasta"), "fasta")]
            df_in = pd.DataFrame(df_data_dict, columns=["id", "seq", "len"]).astype({"len": int})

            id_pattern = r"^(\w+)_split_([\w.]+of\d+):(\d+):(\d+):(\w+):(\w+)$"

            def revise_id(match):
                seqid, segment_position, start, end, TIR_pattern, TSD = match.groups()
                offset = retrieve_split_sequence_offset(segment_position, split_seq_len, overlap_seq_len)
                return f"{seqid}:{int(start) + offset}:{int(end) + offset}:{TIR_pattern}:{TSD}"

            df_in["id"] = df_in["id"].str.replace(id_pattern, revise_id, regex=True)

            df_list.append(df_in.copy())
        except FileNotFoundError:
            continue
        if not flag_debug:
            subprocess.Popen(["rm", "-rf", f])
    return pd.concat(df_list).drop_duplicates(ignore_index=True).sort_values("id", ignore_index=True)


def get_GRF_result_df_native(genome_name, flag_debug):
    GRF_result_dir_name = f"{genome_name}_GRFmite"
    GRF_result_file_path = os.path.join(GRF_result_dir_name, "candidate.fasta")

    # df_data_dict = {"id": [rec.id for rec in SeqIO.parse(GRF_result_file_path, "fasta")],
    # df_data_dict = {"id": list(SeqIO.index(GRF_result_file_path, "fasta")),
    #                 "seq": [str(rec.seq) for rec in SeqIO.parse(GRF_result_file_path, "fasta")],
    #                 "len": [len(rec) for rec in SeqIO.parse(GRF_result_file_path, "fasta")]}
    df_data_dict = [{"id": rec.id, "seq": str(rec.seq), "len": len(rec)}
                    for rec in SeqIO.parse(GRF_result_file_path, "fasta")]

    if not flag_debug:
        subprocess.Popen(["rm", "-rf", GRF_result_dir_name])
    return pd.DataFrame(df_data_dict, columns=["id", "seq", "len"]).astype({"len": int})


def execute(TIRLearner_instance):
    genome_file = TIRLearner_instance.genome_file_path
    genome_name = TIRLearner_instance.genome_name
    TIR_length = TIRLearner_instance.TIR_length
    cpu_cores = TIRLearner_instance.cpu_cores
    flag_debug = TIRLearner_instance.flag_debug
    GRF_path = TIRLearner_instance.GRF_path
    additional_args = TIRLearner_instance.additional_args
    fasta_files_path_list = TIRLearner_instance.split_fasta_files_path_list

    if NO_PARALLEL in additional_args:
        return run_GRF_native(genome_file, genome_name, cpu_cores, TIR_length, flag_debug, GRF_path)
    else:
        return run_GRF_boost(genome_file, genome_name, cpu_cores, TIR_length, flag_debug, GRF_path,
                             fasta_files_path_list)
