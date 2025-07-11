#!/usr/app/env python3
# -*- coding: utf-8 -*-

from shared import *


def TA_repeats_check(df_in: pd.DataFrame, column: str = "seq", percent_threshold: float = 0.7) -> pd.DataFrame:
    df = df_in.copy()
    # df['TA_repeats_perc'] = (df["seq"].str.count('T') + df["seq"].str.count('A')) / df["seq"].str.len()
    df["TA_repeats_check"] = (df[column].str.count('T') + df[column].str.count('A') >=
                              df[column].str.len() * percent_threshold)
    return df[~df["TA_repeats_check"]].drop(columns="TA_repeats_check").reset_index(drop=True)


def _check_N(s: str) -> bool:
    n = s.count("N")
    if n > 0:
        return True
    return False


def _check_N_percentage(s: str) -> bool:
    n = s.count("N")
    p = n / len(s)
    if p >= 0.20:
        return True
    return False


def _find_digits_sum(string: str) -> int:
    pattern = r"(\d+)"
    l = re.findall(pattern, string)
    return sum([int(i) for i in l])


def _TSD_check(x: pd.Series) -> bool:
    TSD = x["TSD"]
    if ((len(TSD) > 6 or 
        TSD == "TAA" or 
        TSD == "TTA" or 
        TSD == "TA" or 
        x["seq"][0:4] == "CACT") or
        x["seq"][0:4] == "GTGA"):
        return True
    return False

#This behavior was pushed to the run_GRF script and this is no longer used at all.
def process_GRF_result(TIRLearner_instance) -> Optional[pd.DataFrame]:
    df_in = TIRLearner_instance.working_df_dict["GRF"]

    if df_in is None:
        print("NOTICE: No TIR candidates was found by GRF.")
        return None

    df = df_in[df_in["len"] >= 50].copy()
    if df.shape[0] == 0:
        print("NOTICE: No valid TIR candidates was found by GRF.")
        return None

    print("  Step 1/7: Getting TIR")
    df["TIR_len"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: _find_digits_sum(x["id"].split(":")[-2]), axis=1)
    df["TIR"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: x["seq"][0:x["TIR_len"]], axis=1)

    print("  Step 2/7: Checking TA repeats on sequence")
    df = TA_repeats_check(df)

    print("  Step 3/7: Checking percentage of N on sequence")
    df["check_N_per_seq_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if _check_N_percentage(x["seq"]) else False, axis=1)
    df = df.dropna(ignore_index=True).copy()

    print("  Step 4/7: Checking TA repeats on TIR")
    df = TA_repeats_check(df, "TIR")

    print("  Step 5/7: Checking N existence on TIR")
    df["check_N_TIR_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if _check_N(x["TIR"]) else False, axis=1)
    df = df.dropna(ignore_index=True).copy()

    print("  Step 6/7: Getting TSD")
    df["TSD"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: x["id"].split(":")[-1], axis=1)

    print("  Step 7/7: Checking TSD")
    df["TSD_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: True if _TSD_check(x) else np.nan, axis=1)
    df = df.dropna(ignore_index=True).loc[:, ["id", "seq"]].copy()

    if df.shape[0] == 0:
        print("[NOTICE] No TIR candidates was found by GRF.")
        return None

    df["id"] = ">" + df["id"]
    return df

#This behavior was pushed to the run_TIRvish and get_fasta_sequence scripts and this is no longer used at all.
def process_TIRvish_result(TIRLearner_instance) -> Optional[pd.DataFrame]:
    df_in = TIRLearner_instance["TIRvish"]

    if df_in is None:
        print("[NOTICE] No TIR candidates was found by TIRvish.")
        return None

    df = df_in[df_in["end"] - df_in["start"] + 1 >= 50].copy()
    if df.shape[0] == 0:
        print("[NOTICE] No valid TIR candidates was found by TIRvish.")
        return None

    print("  Step 1/5: Getting TIR")
    df["TIR1_start"] = df["TIR1_start"] - df["start"]
    df.loc[df["TIR1_start"] < 0, "TIR1_start"] = 0
    df["TIR1_end"] = df["TIR1_end"] - df["start"]

    df["TIR2_start"] = df["TIR2_start"] - df["start"]
    df.loc[df["TIR2_start"] < 0, "TIR2_start"] = 0
    df["TIR2_end"] = df["TIR2_end"] - df["start"]

    df["TIR1"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: x["seq"][x["TIR1_start"]: x["TIR1_end"]], axis=1)
    df["TIR2"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: x["seq"][x["TIR2_start"]: x["TIR2_end"]], axis=1)

    print("  Step 2/5: Checking TA repeats on sequence")
    df = TA_repeats_check(df)

    print("  Step 3/5: Checking percentage of N on sequence")
    df["check_N_per_seq_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if _check_N_percentage(x["seq"]) else False, axis=1)
    df = df.dropna(ignore_index=True).copy()

    print("  Step 4/5: Checking TA repeats on TIR")
    df = TA_repeats_check(df, "TIR1")
    df = TA_repeats_check(df, "TIR2")

    print("  Step 5/5: Checking N existence on TIR")
    df["check_N_TIR1_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if _check_N(x["TIR1"]) else False, axis=1)
    df = df.dropna(ignore_index=True).copy()
    df["check_N_TIR2_check"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: np.nan if _check_N(x["TIR2"]) else False, axis=1)
    
    #Ultimately drop everything else.
    df = df.dropna(ignore_index=True).loc[:, ["id", "seq"]].copy()

    if df.shape[0] == 0:
        print("NOTICE: No TIR was found by TIRvish.")
        return None
    return df

#No longer used
def combine_de_novo_result(TIRLearner_instance):
    try:
        df = pd.concat((TIRLearner_instance.get("TIRvish"), TIRLearner_instance.get("GRF")), ignore_index=True)
    except ValueError:
        df = None

    if df is None or df.shape[0] == 0:
        raise SystemExit("[WARN] No TIR was found by TIRvish and GRF in the de novo process. "
                         "Check your input file as well as arguments.")
    # TODO add and use safe exit function that clears temp dir

    df.to_csv(TIRLearner_instance.processed_de_novo_result_file_name, sep='\n', header=False, index=False)
    del TIRLearner_instance["TIRvish"]
    del TIRLearner_instance["GRF"]
