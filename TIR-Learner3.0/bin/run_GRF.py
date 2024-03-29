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
import os

from const import *

from run_TIRvish import retrieve_split_sequence_offset


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

    # add_three_spaces = r" | sed 's/^/   /'"
    # subprocess.Popen(grf + shell_filter + add_three_spaces, shell=True).wait()
    subprocess.Popen(grf + shell_filter, shell=True).wait()


def GRF_mp(genome_file_path, genome_name, cpu_cores, TIR_length, GRF_path):
    GRF_working_dir = os.path.dirname(genome_file_path)
    genome_file_name = os.path.basename(genome_file_path)
    os.chdir(GRF_working_dir)
    GRF(genome_file_name, genome_name, cpu_cores, TIR_length, GRF_path)
    os.chdir("../")

# def run_GRF_native(filtered_genome_file_name, GRF_path, cpu_cores, TIR_length):
#     GRF(GRF_path, filtered_genome_file_name, int(cpu_cores * thread_core_ratio), TIR_length)
#     subprocess.Popen(["unlink", filtered_genome_file_name])


def cpu_cores_allocation_GRF_boost(cpu_cores, num_seq):
    num_threads = int(cpu_cores * thread_core_ratio)
    num_process = int(math.sqrt(num_threads/16) + 1) * int(1 + num_seq/2)
    # num_process = floor((num_threads/16)^(1/2) + 1) * floor(1 + num_seq/2)
    return num_process, num_threads


# def run_GRF_boost(fasta_files_list, filtered_genome_file_name, num_seq, GRF_path, cpu_cores, TIR_length):
def run_GRF_boost(genome_name, cpu_cores, TIR_length, flag_debug, GRF_path, fasta_files_path_list):
    os.makedirs(f"{splited_fasta_tag}_mp", exist_ok=True)  # TODO revise to include mode detection function
    os.chdir(f"./{splited_fasta_tag}_mp")

    print("  Step 1/2: Executing GRF\n")
    mp_args_list = [(file_path, genome_name, 2, TIR_length, GRF_path) for file_path in fasta_files_path_list]
    # TODO cpu_cores algorithm needed here
    with mp.Pool(cpu_cores) as pool:
        pool.starmap(GRF_mp, mp_args_list)
    print()

    print("  Step 2/2: Getting GRF result")
    return get_GRF_result_df_boost(fasta_files_path_list, genome_name, flag_debug,
                                   TIRvish_split_seq_len, TIRvish_overlap_seq_len)


    # collect_results(fasta_files_path_list, genome_name)
    # subprocess.Popen(["unlink", filtered_genome_file_name])
    # subprocess.Popen(["rm", "-rf", "GRFmite_mp"])


# def collect_results(split_fasta_files_path_list, genome_name):
#     GRF_result_dir_name = f"{genome_name}_GRFmite"
#     GRF_result_dir_list = [os.path.join(os.path.dirname(file), GRF_result_dir_name) for file in
#                            split_fasta_files_path_list]
#     os.makedirs(GRF_result_dir_name, exist_ok=True)
#     with open(os.path.join(GRF_result_dir_name, "candidate.fasta"), "wb") as des:
#         for f in GRF_result_dir_list:
#             try:
#                 with open(os.path.join(f, "candidate.fasta"), "rb") as src:
#                     shutil.copyfileobj(src, des)
#             except FileNotFoundError:
#                 continue


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
    cpu_cores = TIRLearner_instance.cpu_cores
    flag_debug = TIRLearner_instance.flag_debug
    GRF_path = TIRLearner_instance.GRF_path
    TIR_length = TIRLearner_instance.TIR_length
    GRF_mode = TIRLearner_instance.GRF_mode
    split_fasta_files_path_list = TIRLearner_instance.split_fasta_files_path_list


    # records_split_file_name, filtered_genome_file_name, num_seq = prepare_fasta(genome_file, genome_name,
    #                                                                             GRF_mode, drop_seq_len)
    # TODO: More algorithms for mode selection here

    # print(os.listdir(os.getcwd()))
    # for i in os.listdir(os.getcwd()):
    #     shutil.copyfile(os.path.join(os.getcwd(), i), os.path.join(TIRLearner_instance.output_dir, i))

    # run_GRF_native(filtered_genome_file_name, GRF_path, cpu_cores, TIR_length)
    # get_GRF_result_df_native(genome_name, flag_debug)

    # print("  Step 1/2: Executing GRF\n")
    return run_GRF_boost(genome_name, cpu_cores, TIR_length, flag_debug, GRF_path, split_fasta_files_path_list)
