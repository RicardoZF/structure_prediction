#!/sugon_store/pub_data/tools/alphafold3/envs/af3/bin/python3.11

# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Library to run Jackhmmer from Python."""

import os
import tempfile
import subprocess
import time
import pathlib
import string
from typing import Any
from collections.abc import Sequence
import multiprocessing

from absl import logging, flags, app

import dataclasses
from typing import Protocol
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from concurrent import futures
  
_DEFAULT_DB_DIR = pathlib.Path('/sugon_store/pub_data/databases/alphafold3')

# Input and output paths.
_FASTA_PATH = flags.DEFINE_string(
    'fasta_path',
    None,
    'Path to the input FASTA file.',
)

_OUTPUT_DIR = flags.DEFINE_string(
    'output_dir',
    None,
    'Path to a directory where the results will be saved.',
)

# Binary paths.
_JACKHMMER_BINARY_PATH = flags.DEFINE_string(
    'jackhmmer_binary_path',
    '/sugon_store/pub_data/tools/alphafold3/envs/af3/bin/jackhmmer',
    'Path to the Jackhmmer binary.',
)

_REFORMAT_BINARY_PATH = flags.DEFINE_string(
    'reformat_binary_path',
    '/sugon_store/pub_data/tools/alphafold3/envs/af3/bin/reformat.pl',
    'Path to the reformat.pl binary (.sto to .a3m.).',
)

# Database paths.
DB_DIR = flags.DEFINE_multi_string(
    'db_dir',
    (_DEFAULT_DB_DIR.as_posix(),),
    'Path to the directory containing the databases. Can be specified multiple'
    ' times to search multiple directories in order.',
)

_SMALL_BFD_DATABASE_PATH = flags.DEFINE_string(
    'small_bfd_database_path',
    '${DB_DIR}/bfd-first_non_consensus_sequences.fasta',
    'Small BFD database path, used for protein MSA search.',
)
_MGNIFY_DATABASE_PATH = flags.DEFINE_string(
    'mgnify_database_path',
    '${DB_DIR}/mgy_clusters_2022_05.fa',
    'Mgnify database path, used for protein MSA search.',
)
_UNIPROT_CLUSTER_ANNOT_DATABASE_PATH = flags.DEFINE_string(
    'uniprot_cluster_annot_database_path',
    '${DB_DIR}/uniprot_all_2021_04.fa',
    'UniProt database path, used for protein paired MSA search.',
)
_UNIREF90_DATABASE_PATH = flags.DEFINE_string(
    'uniref90_database_path',
    '${DB_DIR}/uniref90_2022_05.fa',
    'UniRef90 database path, used for MSA search. The MSA obtained by '
    'searching it is used to construct the profile for template search.',
)

# Number of CPUs to use for MSA tools.
_JACKHMMER_N_CPU = flags.DEFINE_integer(
    'jackhmmer_n_cpu',
    min(multiprocessing.cpu_count(), 8),
    'Number of CPUs to use for Jackhmmer. Default to min(cpu_count, 8). Going'
    ' beyond 8 CPUs provides very little additional speedup.',
)

def replace_db_dir(path_with_db_dir: str, db_dirs: Sequence[str]) -> str:
  """Replaces the DB_DIR placeholder in a path with the given DB_DIR."""
  template = string.Template(path_with_db_dir)
  if 'DB_DIR' in template.get_identifiers():
    for db_dir in db_dirs:
      path = template.substitute(DB_DIR=db_dir)
      if os.path.exists(path):
        return path
    raise FileNotFoundError(
        f'{path_with_db_dir} with ${{DB_DIR}} not found in any of {db_dirs}.'
    )
  if not os.path.exists(path_with_db_dir):
    raise FileNotFoundError(f'{path_with_db_dir} does not exist.')
  return path_with_db_dir


def _keep_line(line: str, seqnames: Set[str]) -> bool:
    """Function to decide which lines to keep."""
    if not line.strip():
        return True
    if line.strip() == "//":  # End tag
        return True
    if line.startswith("# STOCKHOLM"):  # Start tag
        return True
    if line.startswith("#=GC RF"):  # Reference Annotation Line
        return True
    if line[:4] == "#=GS":  # Description lines - keep if sequence in list.
        _, seqname, _ = line.split(maxsplit=2)
        return seqname in seqnames
    elif line.startswith("#"):  # Other markup - filter out
        return False
    else:  # Alignment data - keep if sequence in list.
        seqname = line.partition(" ")[0]
        return seqname in seqnames

def truncate_stockholm_msa(stockholm_msa_path: str, max_sequences: int) -> str:
    """Reads + truncates a Stockholm file while preventing excessive RAM usage."""
    seqnames = set()
    filtered_lines = []

    with open(stockholm_msa_path) as f:
        for line in f:
            if line.strip() and not line.startswith(("#", "//")):
                # Ignore blank lines, markup and end symbols - remainder are alignment
                # sequence parts.
                seqname = line.partition(" ")[0]
                seqnames.add(seqname)
                if len(seqnames) >= max_sequences:
                    break

        f.seek(0)
        for line in f:
            if _keep_line(line, seqnames):
                filtered_lines.append(line)

    return "".join(filtered_lines)


def _convert_sto_seq_to_a3m(
    query_non_gaps: Sequence[bool], sto_seq: str
) -> Iterable[str]:
    for is_query_res_non_gap, sequence_res in zip(query_non_gaps, sto_seq):
        if is_query_res_non_gap:
            yield sequence_res
        elif sequence_res != "-":
            yield sequence_res.lower()


def convert_stockholm_to_a3m(
    stockholm_format: str,
    max_sequences: Optional[int] = None,
    remove_first_row_gaps: bool = True,
) -> str:
    """Converts MSA in Stockholm format to the A3M format."""
    descriptions = {}
    sequences = {}
    reached_max_sequences = False

    for line in stockholm_format.splitlines():
        reached_max_sequences = max_sequences and len(sequences) >= max_sequences
        if line.strip() and not line.startswith(("#", "//")):
            # Ignore blank lines, markup and end symbols - remainder are alignment
            # sequence parts.
            seqname, aligned_seq = line.split(maxsplit=1)
            if seqname not in sequences:
                if reached_max_sequences:
                    continue
                sequences[seqname] = ""
            sequences[seqname] += aligned_seq

    for line in stockholm_format.splitlines():
        if line[:4] == "#=GS":
            # Description row - example format is:
            # #=GS UniRef90_Q9H5Z4/4-78            DE [subseq from] cDNA: FLJ22755 ...
            columns = line.split(maxsplit=3)
            seqname, feature = columns[1:3]
            value = columns[3] if len(columns) == 4 else ""
            if feature != "DE":
                continue
            if reached_max_sequences and seqname not in sequences:
                continue
            descriptions[seqname] = value
            if len(descriptions) == len(sequences):
                break

    # Convert sto format to a3m line by line
    a3m_sequences = {}
    if remove_first_row_gaps:
        # query_sequence is assumed to be the first sequence
        query_sequence = next(iter(sequences.values()))
        query_non_gaps = [res != "-" for res in query_sequence]
    for seqname, sto_sequence in sequences.items():
        # Dots are optional in a3m format and are commonly removed.
        out_sequence = sto_sequence.replace(".", "")
        if remove_first_row_gaps:
            out_sequence = "".join(
                _convert_sto_seq_to_a3m(query_non_gaps, out_sequence)
            )
        a3m_sequences[seqname] = out_sequence

    fasta_chunks = (
        f">{k} {descriptions.get(k, '')}\n{a3m_sequences[k]}" for k in a3m_sequences
    )
    return "\n".join(fasta_chunks) + "\n"  # Include terminating newline.


class SubprocessUtils():
    @classmethod
    def create_query_fasta_file(cls, sequence: str, path: str, linewidth: int = 80):
        """Creates a fasta file with the sequence with line width limit."""
        with open(path, 'w') as f:
            f.write('>query\n')
            i = 0
            while i < len(sequence):
                f.write(f'{sequence[i:(i + linewidth)]}\n')
                i += linewidth

    @classmethod
    def check_binary_exists(cls, path: str, name: str) -> None:
        """Checks if a binary exists on the given path and raises otherwise."""
        if not os.path.exists(path):
            raise RuntimeError(f'{name} binary not found at {path}')

    @classmethod
    def run(cls, 
        cmd: Sequence[str],
        cmd_name: str,
        log_on_process_error: bool = False,
        log_stderr: bool = False,
        log_stdout: bool = False,
        max_out_streams_len: int | None = 500_000,
        **run_kwargs,
    ) -> subprocess.CompletedProcess[Any]:
        """Launches a subprocess, times it, and checks for errors.

        Args:
            cmd: Command to launch.
            cmd_name: Human-readable command name to be used in logs.
            log_on_process_error: Whether to use `logging.error` to log the process'
            stderr on failure.
            log_stderr: Whether to log the stderr of the command.
            log_stdout: Whether to log the stdout of the command.
            max_out_streams_len: Max length of prefix of stdout and stderr included in
            the exception message. Set to `None` to disable truncation.
            **run_kwargs: Any other kwargs for `subprocess.run`.

        Returns:
            The completed process object.

        Raises:
            RuntimeError: if the process completes with a non-zero return code.
        """

        logging.info('Launching subprocess "%s"', ' '.join(cmd))

        start_time = time.time()
        try:
            completed_process = subprocess.run(
                cmd,
                check=True,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
                text=True,
                **run_kwargs,
            )
        except subprocess.CalledProcessError as e:
            if log_on_process_error:
                # Logs have a 15k character limit, so log the error line by line.
                logging.error('%s failed. %s stderr begin:', cmd_name, cmd_name)
            for error_line in e.stderr.splitlines():
                if stripped_error_line := error_line.strip():
                    logging.error(stripped_error_line)
            logging.error('%s stderr end.', cmd_name)

            error_msg = (
                f'{cmd_name} failed'
                f'\nstdout:\n{e.stdout[:max_out_streams_len]}\n'
                f'\nstderr:\n{e.stderr[:max_out_streams_len]}'
            )
            raise RuntimeError(error_msg) from e
        end_time = time.time()

        logging.info('Finished %s in %.3f seconds', cmd_name, end_time - start_time)
        stdout, stderr = completed_process.stdout, completed_process.stderr

        if log_stdout and stdout:
            logging.info('%s stdout:\n%s', cmd_name, stdout)

        if log_stderr and stderr:
            logging.info('%s stderr:\n%s', cmd_name, stderr)

        return completed_process

@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class MsaToolResult:
    """The result of a MSA tool query."""

    target_sequence: str
    e_value: float
    a3m: str


class MsaTool(Protocol):
    """Interface for MSA tools."""

    def query(self, target_sequence: str) -> MsaToolResult:
        """Runs the MSA tool on the target sequence."""

class Jackhmmer(MsaTool):
    """Python wrapper of the Jackhmmer binary."""

    def __init__(
        self,
        *,
        binary_path: str,
        #database_path: str,
        n_cpu: int = 8,
        n_iter: int = 3,
        e_value: float | None = 1e-3,
        z_value: float | int | None = None,
        max_sequences: int = 5000,
        filter_f1: float = 5e-4,
        filter_f2: float = 5e-5,
        filter_f3: float = 5e-7,
    ):
        """Initializes the Python Jackhmmer wrapper.

        Args:
        binary_path: The path to the jackhmmer executable.
        database_path: The path to the jackhmmer database (FASTA format).
        n_cpu: The number of CPUs to give Jackhmmer.
        n_iter: The number of Jackhmmer iterations.
        e_value: The E-value, see Jackhmmer docs for more details.
        z_value: The Z-value representing the number of comparisons done (i.e
            correct database size) for E-value calculation.
        max_sequences: Maximum number of sequences to return in the MSA.
        filter_f1: MSV and biased composition pre-filter, set to >1.0 to turn off.
        filter_f2: Viterbi pre-filter, set to >1.0 to turn off.
        filter_f3: Forward pre-filter, set to >1.0 to turn off.

        Raises:
        RuntimeError: If Jackhmmer binary not found within the path.
        """
        self.binary_path = binary_path

        SubprocessUtils.check_binary_exists(
            path=self.binary_path, name='Jackhmmer'
        )

        self.n_cpu = n_cpu
        self.n_iter = n_iter
        self.e_value = e_value
        self.z_value = z_value
        self.max_sequences = max_sequences
        self.filter_f1 = filter_f1
        self.filter_f2 = filter_f2
        self.filter_f3 = filter_f3

    def query(self, target_sequence: str, database_path: str, output_a3m: str = None) -> MsaToolResult:
        """Queries the database using Jackhmmer."""
        if not os.path.exists(database_path):
            raise ValueError(f'Could not find Jackhmmer database {database_path}')
        logging.info('Query sequence: %s', target_sequence)
        with tempfile.TemporaryDirectory() as query_tmp_dir:
            input_fasta_path = os.path.join(query_tmp_dir, 'query.fasta')
            SubprocessUtils.create_query_fasta_file(
                sequence=target_sequence, path=input_fasta_path
            )

            output_sto_path = os.path.join(query_tmp_dir, 'output.sto')
            # TODO test
            # output_sto_path = 'test.sto'
            
            # The F1/F2/F3 are the expected proportion to pass each of the filtering
            # stages (which get progressively more expensive), reducing these
            # speeds up the pipeline at the expensive of sensitivity.  They are
            # currently set very low to make querying Mgnify run in a reasonable
            # amount of time.
            cmd_flags = [
                *('-o', '/dev/null'),  # Don't pollute stdout with Jackhmmer output.
                *('-A', output_sto_path),
                '--noali',
                *('--F1', str(self.filter_f1)),
                *('--F2', str(self.filter_f2)),
                *('--F3', str(self.filter_f3)),
                *('--cpu', str(self.n_cpu)),
                *('-N', str(self.n_iter)),
            ]
            
            # Report only sequences with E-values <= x in per-sequence output.
            if self.e_value is not None:
                cmd_flags.extend(['-E', str(self.e_value)])

                # Use the same value as the reporting e-value (`-E` flag).
                cmd_flags.extend(['--incE', str(self.e_value)])

            if self.z_value is not None:
                cmd_flags.extend(['-Z', str(self.z_value)])

            cmd = (
                [self.binary_path]
                + cmd_flags
                + [input_fasta_path, database_path]
            )

            SubprocessUtils.run(
                cmd=cmd,
                cmd_name='Jackhmmer',
                log_stdout=False,
                log_stderr=True,
                log_on_process_error=True,
            )
            
            if self.max_sequences is None:
                with open(output_sto_path) as f:
                    sto = f.read()
            else:
                sto = truncate_stockholm_msa(output_sto_path, max_sequences=self.max_sequences)
            
            # print(sto)
            a3m = convert_stockholm_to_a3m(sto)
            
            output_a3m = os.path.join(query_tmp_dir, 'output.a3m') if not output_a3m else output_a3m
            with open(output_a3m, 'w') as f:
                f.write(output_a3m)
            # sto_to_a3m_cmd = [_REFORMAT_BINARY_PATH.value, output_sto_path, output_a3m]
            # SubprocessUtils.run(
            #     cmd=sto_to_a3m_cmd,
            #     cmd_name='.sto to .a3m',
            #     log_stdout=False,
            #     log_stderr=True,
            #     log_on_process_error=True,
            #     )
            
            # with open(output_a3m) as f:
            #     a3m = f.read()
                
        return MsaToolResult(
            target_sequence=target_sequence, a3m=a3m, e_value=self.e_value
        )

def main(_):
    # Make sure we can create the output directory before running anything.
    try:
        os.makedirs(_OUTPUT_DIR.value, exist_ok=True)
    except OSError as e:
        print(f'Failed to create output directory {_OUTPUT_DIR.value}: {e}')
        raise
    
    jackhmmer_n_cpu = _JACKHMMER_N_CPU.value
    
    # fasta to sequences
    sequences = []
    msa_position = []
    with open(_FASTA_PATH.value) as f:
        ndx = 0
        for line in f:
            if not line.startswith('>'):
                seq = line.strip()
                msa_position.append(ndx)
                if seq not in sequences:
                    sequences.append(seq)
                    ndx += 1
    
    uniref90_msa_jackhmmer = Jackhmmer(
        binary_path = _JACKHMMER_BINARY_PATH.value,
        n_cpu=jackhmmer_n_cpu,
        n_iter=1,
        e_value=1e-4,
        z_value=None,
        max_sequences=10_000,
    )
    
    mgnify_msa_jackhmmer = Jackhmmer(
        binary_path = _JACKHMMER_BINARY_PATH.value,
        n_cpu=jackhmmer_n_cpu,
        n_iter=1,
        e_value=1e-4,
        z_value=None,
        max_sequences=5_000,
    )
    
    small_bfd_msa_jackhmmer = Jackhmmer(
        binary_path = _JACKHMMER_BINARY_PATH.value,
        n_cpu=jackhmmer_n_cpu,
        n_iter=1,
        e_value=1e-4,
        # Set z_value=138_515_945 to match the z_value used in the paper.
        # In practice, this has minimal impact on predicted structures.
        z_value=None,
        max_sequences=5_000,
    )
    
    uniprot_msa_jackhmmer = Jackhmmer(
        binary_path = _JACKHMMER_BINARY_PATH.value,
        n_cpu=jackhmmer_n_cpu,
        n_iter=1,
        e_value=1e-4,
        z_value=None,
        max_sequences=50_000,
    )
    
    msa_io = []
    for sequence in sequences:
        """Processes a single protein chain."""
        logging.info('Getting protein MSAs for sequence %s', sequence)
        msa_start_time = time.time()
        expand_path = lambda x: replace_db_dir(x, DB_DIR.value) 
        # Run various MSA tools in parallel. Use a ThreadPoolExecutor because
        # they're not blocked by the GIL, as they're sub-shelled out.
        with futures.ThreadPoolExecutor(max_workers=4) as executor:
            uniref90_msa_future = executor.submit(
                uniref90_msa_jackhmmer.query,
                target_sequence = sequence,
                database_path = expand_path(_UNIREF90_DATABASE_PATH.value),
            )
            mgnify_msa_future = executor.submit(
                mgnify_msa_jackhmmer.query,
                target_sequence = sequence,
                database_path = expand_path(_MGNIFY_DATABASE_PATH.value),
            )
            small_bfd_msa_future = executor.submit(
                small_bfd_msa_jackhmmer.query,
                target_sequence = sequence,
                database_path = expand_path(_SMALL_BFD_DATABASE_PATH.value),
            )
            uniprot_msa_future = executor.submit(
                uniprot_msa_jackhmmer.query,
                target_sequence = sequence,
                database_path = expand_path(_UNIPROT_CLUSTER_ANNOT_DATABASE_PATH.value),
            )
        logging.info(
            'Getting protein MSAs took %.2f seconds for sequence %s',
            time.time() - msa_start_time,
            sequence,
        )
        
        # test 
        # sequence_data =  small_bfd_msa_future.result().a3m
        sequence_data = "".join({uniref90_msa_future.result().a3m, mgnify_msa_future.result().a3m, small_bfd_msa_future.result().a3m, uniprot_msa_future.result().a3m})
        msa_io.append(sequence_data)
    
    # output .a3m
    for index in msa_position:
        out_msa = os.path.join(_OUTPUT_DIR.value, f"{index}.a3m")
        with open(out_msa, 'w') as f:
            f.write(msa_io[index])

if __name__ == '__main__':
    flags.mark_flags_as_required([
        'fasta_path', 'output_dir'
    ])
    app.run(main)
