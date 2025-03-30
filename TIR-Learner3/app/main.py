#!/usr/app/env python3
# -*- coding: utf-8 -*-
# Tianyu (Sky) Lu (tianyu@lu.fm)
# 2025-03-29

from shared import *

import blast_reference
import process_homology
import run_TIRvish
import run_GRF
import process_de_novo_result
import prepare_data
import CNN_predict
import get_fasta_sequence
import check_TIR_TSD
import post_processing


def get_timestamp_now_utc_iso8601(date_sep: str = '-', time_sep: str = ':') -> str:
    """
    Generate current UTC timestamp in ISO 8601 format with customizable separators.

    The timestamp follows the ISO 8601 standard (YYYY-MM-DDTHH:MM:SSZ) but allows
    customization of the date and time component separators.

    :param date_sep: Separator character between date components (year, month, day).
                    Defaults to '-' (ISO 8601 standard).
    :type date_sep: str
    :param time_sep: Separator character between time components (hour, minute, second).
                    Defaults to ':' (ISO 8601 standard).
    :type time_sep: str
    :return: Formatted UTC timestamp string in ISO 8601 format with specified separators
    :rtype: str
    :raises: Does not raise any exceptions

    Example::

        >>> get_timestamp_now_utc_iso8601()
        "2023-10-05T14:30:00Z"

        >>> get_timestamp_now_utc_iso8601(time_sep='-')
        "2023-10-05T14-30-00Z"

        >>> get_timestamp_now_utc_iso8601(date_sep='', time_sep='')
        "20231005T143000Z"
    """
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        f"%Y{date_sep}%m{date_sep}%dT%H{time_sep}%M{time_sep}%SZ"
    )


def humanize_float(number: float, base: Union[float, int], sig_fig: int, units: Sequence[str]) -> str:
    """
    Convert a floating-point number to human-readable format with scaling.

    Scales the input number according to the specified base and returns a string
    representation with the appropriate unit and significant figures.

    :param number: Floating-point number to convert
    :type number: float
    :param base: Base value for scaling (e.g., 1000 for SI units)
    :type base: Union[float, int]
    :param sig_fig: Number of significant figures to include in output
    :type sig_fig: int
    :param units: Sequence of unit labels for different scales (e.g., ("B", "KB", "MB"))
    :type units: Sequence[str]
    :return: Human-readable string representation of the number
    :rtype: str

    Example::

        >>> humanize_float(1500, 1000, 3, ("B", "KB", "MB"))
        "1.50 KB"
    """
    return next(f"{number / base ** i:.{sig_fig}g} {unit}" for
                i, unit in enumerate(units) if number < base ** (i + 1) or i == len(units) - 1)


def humanize_time(seconds: float) -> str:
    """
    Convert seconds to a human-readable time format.

    Scales the input seconds into appropriate units (seconds, minutes, hours)
    and formats with proper significant figures.

    :param seconds: Time duration in seconds to convert
    :type seconds: float
    :return: Human-readable time string
    :rtype: str

    Example::

        >>> humanize_time(125)
        "2.1 min"
        >>> humanize_time(3600)
        "1.0 h"
    """
    return humanize_float(seconds, 60, 2 + PRINT_DECIMAL_SIGNIFICANT_FIGURES, ("s", "min", "h"))


def humanize_file_size(file_size: float) -> str:
    """
    Convert file size in bytes to a human-readable format.

    Scales the input size through binary prefixes (KiB, MiB, GiB, etc.)
    and formats with proper significant figures.

    :param file_size: File size in bytes to convert
    :type file_size: float
    :return: Human-readable file size string
    :rtype: str

    Example::

        >>> humanize_file_size(2048)
        "2.0 KiB"
        >>> humanize_file_size(1500000)
        "1.4 MiB"
    """
    return humanize_float(file_size, 1024, 4 + PRINT_DECIMAL_SIGNIFICANT_FIGURES,
                          ("Byte", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB"))


def humanize_sequence_length(length: float) -> str:
    """
    Convert sequence length in base pairs to a human-readable format.

    Scales the input length through metric prefixes (kbp, Mbp, Gbp)
    and formats with proper significant figures.

    :param length: Sequence length in base pairs to convert
    :type length: float
    :return: Human-readable sequence length string
    :rtype: str

    Example::

        >>> humanize_sequence_length(1500)
        "1.5 kbp"
        >>> humanize_sequence_length(2500000)
        "2.5 Mbp"
    """
    return humanize_float(length, 1000, 4 + PRINT_DECIMAL_SIGNIFICANT_FIGURES, ("bp", "kbp", "Mbp", "Gbp"))


def get_process_memory_usage_in_byte(pid: int = os.getpid(),
                                     detail: bool = True) -> Optional[dict[str, Union[dict[str, int], int]]]:
    """
    Retrieve memory usage for a specific process and all its child processes.

    :param pid: Process ID to analyze. Defaults to the current process ID.
    :type pid: int
    :param detail: If True, returns a detailed dictionary for each process.
                   If False, returns only the total memory usage.
    :type detail: bool
    :return: A dictionary containing memory usage information.
             If `detail` is True, the dictionary contains individual process details and a total.
             If `detail` is False, the dictionary contains only the total memory usage.
    :rtype: Optional[dict[str, Union[dict[str, int], int]]]
    """

    try:
        parent = psutil.Process(pid)
        result: Dict[Union[int, str], Dict[str, int]] = {}

        mem = parent.memory_info()
        total_rss: int = mem.rss
        total_vms: int = mem.vms

        if detail:
            result[pid] = {"rss": mem.rss, "vms": mem.vms}

        try:
            children = parent.children(recursive=True)

            for child in children:
                try:
                    mem = child.memory_info()
                    if detail:
                        result[child.pid] = {"rss": mem.rss, "vms": mem.vms}
                    total_rss += mem.rss
                    total_vms += mem.vms
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass

        if detail:
            result["total"] = {"rss": total_rss, "vms": total_vms}
            return result
        return {"rss": total_rss, "vms": total_vms}

    except (psutil.NoSuchProcess, psutil.AccessDenied) as e:
        print(f"Error accessing process: {e}")
        return None


def get_df_total_memory_usage_in_byte(obj: Optional[Union[dict[str, pd.DataFrame], pd.DataFrame]]) -> int:
    """Calculate the total memory usage of pandas DataFrame(s) in bytes.

    Recursively calculates memory usage for either a single DataFrame or a dictionary of DataFrames.
    Handles None input by returning 0 bytes usage.

    :param obj: Input object containing DataFrame(s). Can be one of:
        - ``None`` (returns 0)
        - Single ``pd.DataFrame`` object
        - Dictionary with string keys and DataFrame values
    :type obj: Optional[Union[dict[str, pd.DataFrame], pd.DataFrame]]

    :return: Total memory usage in bytes. Returns 0 for None input.
    :rtype: int

    :raises ValueError: If input type is not None, DataFrame, or dict of DataFrames

    Example:
        >>> get_df_total_memory_usage_in_byte(pd.DataFrame({'a': [1,2]}))
        56
        >>> get_df_total_memory_usage_in_byte({'df1': pd.DataFrame({'a': [1,2]}), 'df2': pd.DataFrame({'b': [3,4]})})
        [sum of both DataFrames' memory usage]
    """
    if obj is None:
        return 0
    elif type(obj) is pd.DataFrame:
        return obj.memory_usage(index=False, deep=True).sum()
    elif type(obj) is dict:
        total_size = 0
        for k, v in obj.items():
            total_size += get_df_total_memory_usage_in_byte(v)
        return total_size
    else:
        raise ValueError(f"[ERROR] Not supported type \"{type(obj)}\"")


class TIRLearner:
    def __init__(self, genome_file_path: str, genome_name: str, species: str, TIR_length: int,
                 processors: int, para_mode: str,
                 working_dir_path: str, output_dir_path: str, checkpoint_dir_input_path: str,
                 flag_verbose: bool, flag_debug: bool, GRF_path: str, gt_path: str,
                 additional_args: Optional[Tuple[str, ...]]):
        self.genome_file_path: str = genome_file_path
        self.genome_name: str = genome_name
        self.species: str = species

        self.TIR_length: int = TIR_length
        self.processors: int = processors
        self.para_mode: str = para_mode

        self.working_dir_path: str = working_dir_path
        self.output_dir_path: str = output_dir_path
        self.checkpoint_dir_input_path: str = checkpoint_dir_input_path

        self.flag_verbose: bool = flag_verbose
        self.flag_debug: bool = flag_debug

        self.GRF_path: str = GRF_path
        self.gt_path: str = gt_path
        self.additional_args: Tuple[str, ...] = additional_args or tuple()

        self.processed_de_novo_result_file_name: str = PROCESSED_DE_NOVO_RESULT_FILE_NAME_FORMAT_STR.format(
            self.genome_name)

        if CHECKPOINT_OFF not in self.additional_args:
            self.checkpoint_dir_output_path: str = os.path.join(
                self.output_dir_path, CHECKPOINT_DIR_NAME_PREFIX + get_timestamp_now_utc_iso8601(time_sep='-'))
            os.makedirs(self.checkpoint_dir_output_path)

        self.genome_file_stat: Dict[str, Union[float, int]] = {"file_size": -0.1, "num": -1,
                                                               "short_seq_num": 0, "short_seq_perc": -0.1,
                                                               "total_len": 0, "avg_len": -1}
        self.current_progress: List[int] = [0, 0]
        self.working_df_dict: Dict[str, Optional[pd.DataFrame]] = {}
        self.split_fasta_files_path_list: List[str] = []

    def __getitem__(self, key: str) -> Optional[pd.DataFrame]:
        return self.working_df_dict[key]

    def __setitem__(self, key: str, value: Optional[pd.DataFrame]):
        self.working_df_dict[key] = value

    def __delitem__(self, key: str) -> int:
        try:
            del self.working_df_dict[key]
            gc.collect()
            return 0
        except KeyError:
            return -1

    def keys(self):
        return self.working_df_dict.keys()

    def items(self):
        return self.working_df_dict.items()

    def get(self, key: str, value=None) -> Optional[pd.DataFrame]:
        return self.working_df_dict.get(key, value)

    def clear(self):
        self.working_df_dict.clear()
        gc.collect()

    # ------------------------------------------------------------------------------------------------------------------

    def execute(self):
        temp_dir: str = self.__mount_working_dir()
        self.__load_checkpoint_file()
        self.__pre_scan_fasta_file()

        # TODO save scanned fasta file in checkpoint
        # TODO revise execute() to put load_checkpoint_file() and pre_scan_fasta_file() in init
        # TODO add cleanup() method to handle the error during init

        if self.species in REFLIB_AVAILABLE_SPECIES:
            self.execute_m1()
            self.execute_m2()
            self.execute_m3()
            raw_result_df_list: List[pd.DataFrame] = [self["m1"], self["m2"], self["m3"]]
        else:
            self.execute_m4()
            raw_result_df_list: List[pd.DataFrame] = [self["m4"]]

        post_processing.execute(self, raw_result_df_list)
        if CHECKPOINT_OFF not in self.additional_args and not self.flag_debug:
            shutil.rmtree(self.checkpoint_dir_output_path)

        # subprocess.Popen(["unlink", self.genome_file_path]).wait()
        # os.rmdir(self.working_dir_path)
        if not self.flag_debug:
            if temp_dir is not None:
                shutil.rmtree(temp_dir)
            else:
                shutil.rmtree(self.working_dir_path)

    def __pre_scan_fasta_file(self):
        print("Doing pre-scan for genome file...")
        start_time: float = time.perf_counter()
        self.genome_file_stat["file_size"] = os.path.getsize(self.genome_file_path)

        try:
            self.genome_file_stat["num"] = len(list(SeqIO.index(self.genome_file_path, "fasta")))
        except Exception:
            raise SystemExit("[ERROR] Duplicate sequence name occurs in the genome file. "
                             "Revise the genome file and try again.")

        records: List[SeqRecord] = []
        for record in SeqIO.parse(self.genome_file_path, "fasta"):
            record.seq = record.seq.upper()
            sequence_str: str = str(record.seq)
            seq_len: int = len(sequence_str)
            drop_seq_len: int = self.TIR_length + 500
            if seq_len < drop_seq_len:
                continue

            nucleotide_diff_set = set(sequence_str) - ACCEPTED_NUCLEOTIDES
            if nucleotide_diff_set:
                print(f"[WARN] Unprocessable nucleotide(s) {nucleotide_diff_set} in sequence {record.id} "
                      f"will be ignored and replaced by \'N\'.")
                record.seq = Seq(re.sub(f"[^{''.join(ACCEPTED_NUCLEOTIDES)}]", "N", sequence_str))

            if SPLITER in record.id:
                print((f"[WARN] Sequence name \"{record.id}\" has reserved string \"{SPLITER}\", "
                       "which makes it incompatible with TIR-Learner and will be replaced with \'_\'."))
                record.id = record.id.replace(SPLITER, "_")

            if len(sequence_str) < SHORT_SEQ_LEN:
                self.genome_file_stat["short_seq_num"] += 1

            self.genome_file_stat["total_len"] += seq_len
            records.append(record)
        checked_genome_file: str = f"{self.genome_name}{SPLITER}checked.fa"
        SeqIO.write(records, checked_genome_file, "fasta")
        self.genome_file_stat["short_seq_perc"] = (self.genome_file_stat["short_seq_num"] /
                                                   self.genome_file_stat["num"])
        self.genome_file_stat["avg_len"] = self.genome_file_stat["total_len"] // self.genome_file_stat["num"]

        end_time: float = time.perf_counter()
        print(f"Genome file scan finished! Time elapsed: {humanize_time(end_time - start_time)}.")
        print(f"  File name: {os.path.basename(self.genome_file_path)}")
        print(f"  File size: {humanize_file_size(self.genome_file_stat['file_size'])}")
        print(f"  Number of sequences: {self.genome_file_stat['num']}")
        print(f"  Number of short sequences: {self.genome_file_stat['short_seq_num']}")
        print(f"  Percentage of short sequences: {self.genome_file_stat['short_seq_perc'] * 100} %")
        print(f"  Average sequence length: {humanize_sequence_length(self.genome_file_stat['avg_len'])}")
        self.genome_file_path: str = os.path.abspath(checked_genome_file)

    def __mount_working_dir(self):
        if self.working_dir_path is None:
            temp_dir: Optional[str] = tempfile.mkdtemp()
            self.working_dir_path = temp_dir
        else:
            temp_dir: Optional[str] = None
            os.makedirs(self.working_dir_path, exist_ok=True)
        self.working_dir_path = os.path.join(self.working_dir_path, SANDBOX_DIR_NAME)
        os.makedirs(self.working_dir_path, exist_ok=True)
        self.working_dir_path = os.path.abspath(self.working_dir_path)
        os.chdir(self.working_dir_path)
        return temp_dir

    # def load_genome_file(self):
    #     genome_file_soft_link = os.path.join(self.execution_dir, "genome_file_soft_link.fa.lnk")
    #     os.symlink(self.genome_file, genome_file_soft_link)
    #     self.genome_file = genome_file_soft_link

    def __get_latest_valid_checkpoint_dir(self, search_dir: str) -> str:
        checkpoint_dirs: List[str] = sorted([f for f in os.listdir(search_dir) if
                                             f.startswith("TIR-Learner_v3_checkpoint_")])
        try:
            checkpoint_dirs.remove(os.path.basename(self.checkpoint_dir_output_path))
        except ValueError:
            pass
        if len(checkpoint_dirs) == 0:
            return "not_found"

        latest_checkpoint_dir_path: str = os.path.join(self.output_dir_path, checkpoint_dirs.pop())
        while not os.listdir(latest_checkpoint_dir_path) and len(checkpoint_dirs) != 0:
            os.rmdir(latest_checkpoint_dir_path)
            latest_checkpoint_dir_path = os.path.join(self.output_dir_path, checkpoint_dirs.pop())

        if not os.listdir(latest_checkpoint_dir_path):
            os.rmdir(latest_checkpoint_dir_path)
            return "not_found"
        return latest_checkpoint_dir_path

    def __load_checkpoint_file(self):
        if CHECKPOINT_OFF in self.additional_args or self.checkpoint_dir_input_path is None:
            return

        if self.checkpoint_dir_input_path == "auto":
            # Search both the output directory and the genome file directory
            self.checkpoint_dir_input_path = self.__get_latest_valid_checkpoint_dir(self.output_dir_path)
            if self.checkpoint_dir_input_path == "not_found":
                genome_file_directory = os.path.dirname(self.genome_file_path)
                self.checkpoint_dir_input_path = self.__get_latest_valid_checkpoint_dir(genome_file_directory)
            if self.checkpoint_dir_input_path == "not_found":
                self.__reset_checkpoint_load_state("Unable to find checkpoint file. ")
                return

        checkpoint_info_file_path: str = os.path.join(self.checkpoint_dir_input_path, "info.txt")

        try:
            with open(checkpoint_info_file_path) as checkpoint_info_file:
                if os.path.getsize(checkpoint_info_file_path) == 0:
                    raise EOFError(checkpoint_info_file_path)
                timestamp_iso8601: str = checkpoint_info_file.readline().rstrip()

                execution_progress_info: List[str] = json.loads(checkpoint_info_file.readline().rstrip())
                species: str = execution_progress_info[0]
                module: int = int(execution_progress_info[1])
                step: int = int(execution_progress_info[2])

                if species != self.species:
                    self.__reset_checkpoint_load_state(f"Species \"{species}\" in checkpoint "
                                                     f"mismatch with species \"{self.species}\" in argument. ")
                    return

                self.current_progress = [module, step]

                working_df_filename_dict: Dict[str, Optional[str]] = json.loads(
                    checkpoint_info_file.readline().rstrip())
                for k, v in working_df_filename_dict.items():
                    if v is None:
                        self[k] = None
                        continue
                    df_file_name: str = os.path.join(self.checkpoint_dir_input_path, v)
                    df_dtype_file_name: str = f"{df_file_name}_dtypes.txt"
                    with open(df_dtype_file_name, 'r') as df_dtype_file:
                        if os.path.getsize(df_dtype_file_name) == 0:
                            raise EOFError(df_dtype_file_name)
                        df_dtypes_dict: Dict[str, str] = {k: eval(v) for k, v in
                                                          json.loads(df_dtype_file.readline().rstrip()).items()}
                    self[k] = pd.read_csv(f"{df_file_name}.csv", sep='\t', header=0,
                                          dtype=df_dtypes_dict, engine='c', memory_map=True)

            processed_de_novo_result_checkpoint_file: str = os.path.join(self.checkpoint_dir_input_path,
                                                                         self.processed_de_novo_result_file_name)
            if os.path.exists(processed_de_novo_result_checkpoint_file):
                shutil.copy(processed_de_novo_result_checkpoint_file,
                            os.path.join(self.working_dir_path, self.processed_de_novo_result_file_name))
                shutil.copy(processed_de_novo_result_checkpoint_file,
                            os.path.join(self.checkpoint_dir_output_path, self.processed_de_novo_result_file_name))

            print("Successfully loaded checkpoint:\n"
                  f"  Time: {timestamp_iso8601}\n"
                  f"  Module: {module}\n"
                  f"  Step: {step}")
        except FileNotFoundError as e:
            self.__reset_checkpoint_load_state("Checkpoint file invalid, "
                                             f"\"{os.path.basename(e.filename)}\" is missing. ")
            return
        except EOFError as e:
            self.__reset_checkpoint_load_state(f"Checkpoint file invalid, \"{os.path.basename(str(e))}\" is empty. ")
            return
        except ValueError:
            self.__reset_checkpoint_load_state(f"Checkpoint file invalid, \"{df_file_name}.csv\" is empty. ")
            return

    def __reset_checkpoint_load_state(self, warn_info: str):
        print("[WARN] " + warn_info + "Will skip loading checkpoint and start from the very beginning.")
        self.current_progress = [0, 0]
        self.clear()

    def __save_checkpoint_file(self):
        if CHECKPOINT_OFF in self.additional_args:
            return

        module: int = self.current_progress[0]
        step: int = self.current_progress[1]
        # checkpoint_file_name = f"module_{module}_step_{step}_{timestamp_now_iso8601}.csv"
        working_df_filename_dict: Dict[str, Optional[str]] = {
            k: f"{k}_module_{module}_step_{step}_{get_timestamp_now_utc_iso8601()}" for k in self.keys()}

        for k, v in self.items():
            if v is None:
                working_df_filename_dict[k] = None
                continue
            df_file_name: str = os.path.join(self.checkpoint_dir_output_path, working_df_filename_dict[k])
            v.to_csv(f"{df_file_name}.csv", index=False, header=True, sep='\t')
            with open(f"{df_file_name}_dtypes.txt", 'w') as f:
                try:
                    f.write(json.dumps(v.loc[0, :].apply(type).apply(str).str[8:-2].
                                       str.replace("numpy", "np").to_dict()) + '\n')
                except KeyError:
                    pass

        with open(os.path.join(self.checkpoint_dir_output_path, "info.txt"), 'w') as f:
            f.write(get_timestamp_now_utc_iso8601() + '\n')
            f.write(json.dumps([self.species, module, step]) + '\n')
            f.write(json.dumps(working_df_filename_dict))
            f.write('\n')

        if not self.flag_debug:
            remove_file_set: set = (set(os.listdir(self.checkpoint_dir_output_path)) -
                                    set(map("{}.csv".format, working_df_filename_dict.values())) -
                                    set(map("{}_dtypes.txt".format, working_df_filename_dict.values())) -
                                    {"info.txt", self.processed_de_novo_result_file_name})
            for f in remove_file_set:
                subprocess.Popen(["unlink", os.path.join(self.checkpoint_dir_output_path, f)],
                                 stderr=subprocess.DEVNULL)

    def __save_processed_de_novo_result_checkpoint_file(self):
        if CHECKPOINT_OFF in self.additional_args:
            return

        shutil.copy(self.processed_de_novo_result_file_name,
                    os.path.join(self.checkpoint_dir_output_path, self.processed_de_novo_result_file_name))
        with open(os.path.join(self.checkpoint_dir_output_path, "info.txt"), 'r') as f:
            lines: List[str] = f.readlines()

        module: int = self.current_progress[0]
        step: int = self.current_progress[1]
        lines[1] = json.dumps((self.species, module, step)) + '\n'

        with open(os.path.join(self.checkpoint_dir_output_path, "info.txt"), 'w') as f:
            f.writelines(lines)

    def __progress_check(self, progress_or_module: Union[list[int], int], step: Optional[int] = None) -> bool:
        if type(progress_or_module) is list:
            executing_module: int = progress_or_module[0]
            executing_step: int = progress_or_module[1]
        else:
            executing_module: int = progress_or_module
            executing_step: int = step

        return (self.checkpoint_dir_input_path is None or
                self.current_progress[0] < executing_module or
                (self.current_progress[0] == executing_module and self.current_progress[1] < executing_step))

    def show_current_memory_usage(self):
        process_memory_usage: Optional[dict[str, Union[dict[str, int], int]]] = get_process_memory_usage_in_byte(
            detail=False)
        if process_memory_usage is None:
            return
        print(f"[Current Memory Usage] "
              f"RSS: {humanize_file_size(process_memory_usage['rss'])}, "
              # f"VMS: {humanize_file_size(process_memory_usage['vms'])}, "
              f"DataFrame: {humanize_file_size(get_df_total_memory_usage_in_byte(self.working_df_dict))}")

    def execute_m1(self):
        terminal_print("Module 1 Begin")
        module: str = "Module1"

        # Module 1, Step 1: Blast reference library in genome file
        current_progress: List[int] = [1, 1]
        if self.__progress_check(current_progress):
            blast_reference.blast_genome_file(self)
            self.current_progress = current_progress
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 1, Step 2: Select 100% coverage entries from blast results
        current_progress = [1, 2]
        if self.__progress_check(current_progress):
            self["base"] = process_homology.select_full_coverage(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 1, Step 3: Making blastDB and get candidate FASTA sequences
        current_progress = [1, 3]
        if self.__progress_check(current_progress):
            print("Module 1, Step 3: Making blastDB and get candidate FASTA sequences")
            self["base"] = get_fasta_sequence.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 1, Step 4: Check TIR and TSD
        current_progress = [1, 4]
        if self.__progress_check(current_progress):
            print("Module 1, Step 4: Check TIR and TSD")
            self["base"] = check_TIR_TSD.execute(self, module)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 1, Step 5: Save module result
        current_progress = [1, 5]
        if self.__progress_check(current_progress):
            self["m1"] = self["base"]
            del self["base"]
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        terminal_print("Module 1 Finished")

    def execute_m2(self):
        terminal_print("Module 2 Begin")
        module: str = "Module2"

        # Module 2, Step 1: Run TIRvish to find inverted repeats
        current_progress: List[int] = [2, 1]
        if self.__progress_check(current_progress) and SKIP_TIRVISH not in self.additional_args:
            print("Module 2, Step 1: Run TIRvish to find inverted repeats")
            self["TIRvish"] = run_TIRvish.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 2, Step 2: Process TIRvish results
        current_progress = [2, 2]
        if self.__progress_check(current_progress) and SKIP_TIRVISH not in self.additional_args:
            print("Module 2, Step 2: Process TIRvish results")
            self["TIRvish"] = process_de_novo_result.process_TIRvish_result(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 2, Step 3: Run GRF to find inverted repeats
        current_progress = [2, 3]
        if self.__progress_check(current_progress) and SKIP_GRF not in self.additional_args:
            print("Module 2, Step 3: Run GRF to find inverted repeats")
            self["GRF"] = run_GRF.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 2, Step 4: Process GRF results
        current_progress = [2, 4]
        if self.__progress_check(current_progress) and SKIP_GRF not in self.additional_args:
            print("Module 2, Step 4: Process GRF results")
            self["GRF"] = process_de_novo_result.process_GRF_result(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 2, Step 5: Combine TIRvish and GRF results
        current_progress = [2, 5]
        if self.__progress_check(current_progress):
            print("Module 2, Step 5: Combine TIRvish and GRF results")
            process_de_novo_result.combine_de_novo_result(self)
            self.current_progress = current_progress
            self.__save_processed_de_novo_result_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 2, Step 6: Blast GRF and TIRvish result in reference library
        current_progress = [2, 6]
        if self.__progress_check(current_progress):
            # Checkpoint saving for this step is currently not available
            blast_reference.blast_de_novo_result(self)
            self.current_progress = current_progress
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 2, Step 7: Select 80% similar entries from blast results
        current_progress = [2, 7]
        if self.__progress_check(current_progress):
            self["base"] = process_homology.select_eighty_similarity(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 2, Step 8: Get FASTA sequences from 80% similarity
        current_progress = [2, 8]
        if self.__progress_check(current_progress):
            print("Module 2, Step 8: Get FASTA sequences from 80% similarity")
            self["base"] = get_fasta_sequence.execute(self)
            self["m2_homo"] = self["base"].copy()
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 2, Step 9: Check TIR and TSD
        current_progress = [2, 9]
        if self.__progress_check(current_progress):
            print("Module 2, Step 9: Check TIR and TSD")
            self["base"] = check_TIR_TSD.execute(self, module)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 2, Step 10: Save module result
        current_progress = [2, 10]
        if self.__progress_check(current_progress):
            print("Module 2, Step 10: Save module result")
            self["m2"] = self["base"]
            del self["base"]
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        terminal_print("Module 2 Finished")

    def execute_m3(self):
        terminal_print("Module 3 Begin")
        module: str = "Module3"

        # Module 3, Step 1: Prepare data
        current_progress: List[int] = [3, 1]
        if self.__progress_check(current_progress):
            print("Module 3, Step 1: Prepare data")
            self["base"] = prepare_data.execute(self, self["m2_homo"])
            del self["m2_homo"]
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 3, Step 2: CNN prediction
        current_progress = [3, 2]
        if self.__progress_check(current_progress):
            print("Module 3, Step 2: CNN prediction")
            self["base"] = CNN_predict.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 3, Step 3: Get FASTA sequences from CNN prediction
        current_progress = [3, 3]
        if self.__progress_check(current_progress):
            print("Module 3, Step 3: Get FASTA sequences from CNN prediction")
            self["base"] = get_fasta_sequence.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 3, Step 4: Check TIR and TSD
        current_progress = [3, 4]
        if self.__progress_check(current_progress):
            print("Module 3, Step 4: Check TIR and TSD")
            self["base"] = check_TIR_TSD.execute(self, module)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 3, Step 5: Save module result
        current_progress = [3, 5]
        if self.__progress_check(current_progress):
            print("Module 3, Step 5: Save module result")
            self["m3"] = self["base"]
            del self["base"]
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        terminal_print("Module 3 Finished")

    def execute_m4(self):
        terminal_print("Module 4 Begin")
        module: str = "Module4"

        # Module 4, Step 1: Run TIRvish to find inverted repeats
        current_progress: List[int] = [4, 1]
        if self.__progress_check(current_progress) and SKIP_TIRVISH not in self.additional_args:
            print("Module 4, Step 1: Run TIRvish to find inverted repeats")
            self["TIRvish"] = run_TIRvish.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 4, Step 2: Process TIRvish results
        current_progress = [4, 2]
        if self.__progress_check(current_progress) and SKIP_TIRVISH not in self.additional_args:
            print("Module 4, Step 2: Process TIRvish results")
            self["TIRvish"] = process_de_novo_result.process_TIRvish_result(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 4, Step 3: Run GRF to find inverted repeats
        current_progress = [4, 3]
        if self.__progress_check(current_progress) and SKIP_GRF not in self.additional_args:
            print("Module 4, Step 3: Run GRF to find inverted repeats")
            self["GRF"] = run_GRF.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 4, Step 4: Process GRF results
        current_progress = [4, 4]
        if self.__progress_check(current_progress) and SKIP_GRF not in self.additional_args:
            print("Module 4, Step 4: Process GRF results")
            self["GRF"] = process_de_novo_result.process_GRF_result(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 4, Step 5: Combine TIRvish and GRF results
        current_progress = [4, 5]
        if self.__progress_check(current_progress):
            print("Module 4, Step 5: Combine TIRvish and GRF results")
            process_de_novo_result.combine_de_novo_result(self)
            self.current_progress = current_progress
            self.__save_processed_de_novo_result_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 4, Step 6: Prepare data
        current_progress = [4, 6]
        if self.__progress_check(current_progress):
            print("Module 4, Step 6: Prepare data")
            self["base"] = prepare_data.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 4, Step 7: CNN prediction
        print("Module 4, Step 7: CNN prediction")
        current_progress = [4, 7]
        if self.__progress_check(current_progress):
            self["base"] = CNN_predict.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 4, Step 8: Get FASTA sequences from CNN prediction
        current_progress = [4, 8]
        if self.__progress_check(current_progress):
            print("Module 4, Step 8: Get FASTA sequences from CNN prediction")
            self["base"] = get_fasta_sequence.execute(self)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 4, Step 9: Check TIR and TSD
        current_progress = [4, 9]
        if self.__progress_check(current_progress):
            print("Module 4, Step 9: Check TIR and TSD")
            self["base"] = check_TIR_TSD.execute(self, module)
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        # Module 4, Step 10: Save module result
        current_progress = [4, 10]
        if self.__progress_check(current_progress):
            self["m4"] = self["base"]
            del self["base"]
            self.current_progress = current_progress
            self.__save_checkpoint_file()
            if self.flag_debug:
                self.show_current_memory_usage()

        terminal_print("Module 4 Finished")
