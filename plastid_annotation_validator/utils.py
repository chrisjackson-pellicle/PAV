#!/usr/bin/env python

import os
import sys
import io
import pstats
import logging
import datetime
import shutil
from Bio.Seq import Seq
import subprocess
from Bio import SeqIO
import re
import threading
from multiprocessing import Manager
import time
import gzip



class LogManager:
    """Global log manager that handles automatic cleanup."""
    
    def __init__(self):
        self.logger = None
        self.log_queue = None
        self.log_listener = None
    
    def setup(self, name, log_file, log_directory=None, **kwargs):
        """Set up the global logger."""
        self.logger, self.log_queue, self.log_listener = setup_queue_logger(
            name, log_file, log_directory, **kwargs
        )
        return self.logger, self.log_queue, self.log_listener
    
    def cleanup(self, timeout=2.0):
        """Clean up the global logger."""
        if self.log_listener:
            self.log_listener.flush(timeout=timeout)
            self.log_listener.stop()
            self.log_listener.join()
    
    def handle_error(self, error, tb, operation_name, id=None):
        """Handle errors with automatic cleanup."""
        if self.logger:
            if id:
                self.logger.error(f'{"[ERROR]:":15} {operation_name} failed for ID {id}')
            else:
                self.logger.error(f'{"[ERROR]:":15} {operation_name} failed')

            if isinstance(error, subprocess.CalledProcessError):
                self.logger.error(f'{"":15} Command: {" ".join(error.cmd)}')
                self.logger.error(f'{"":15} Return code: {error.returncode}')
                self.logger.error(f'{"":15} STDOUT: {error.stdout}')
                self.logger.error(f'{"":15} STDERR: {error.stderr}')
                self.logger.error(f'{"":15} Full traceback:')
                self.logger.error(tb)
            
            elif isinstance(error, Exception):
                self.logger.error(f'{"":15} Exception: {error}')
                self.logger.error(f'{"":15} Full traceback:')
                self.logger.error(tb)

        self.cleanup()
        sys.exit(1)

# Global instance
log_manager = LogManager()

class QueueHandler(logging.Handler):
    """A handler that sends log records to a queue for multiprocessing-safe logging.

    This handler is designed to work with multiprocessing environments where
    direct logging to files or console from worker processes can cause issues.
    Instead, it sends log records to a shared queue that is consumed by a
    listener thread in the main process.

    Attributes:
        queue: A multiprocessing-safe queue for sending log records
    """

    def __init__(self, queue):
        """Initialize the QueueHandler.

        Args:
            queue: A multiprocessing-safe queue (e.g., from multiprocessing.Manager)
        """
        super().__init__()
        self.queue = queue

    def emit(self, record):
        """Send the log record to the queue.

        Args:
            record: The log record to send
        """
        try:
            self.queue.put_nowait(record)
        except Exception:
            self.handleError(record)


class QueueListener(threading.Thread):
    """A thread that listens for log records from a queue and forwards them to handlers.

    This listener runs in the main process and consumes log records from a shared
    queue, forwarding them to appropriate handlers (file, console, etc.) based on
    their log levels. This ensures thread-safe logging in multiprocessing environments.

    Attributes:
        queue: A multiprocessing-safe queue for receiving log records
        handlers: List of logging handlers to forward records to
        _stop_event: Threading event to signal when to stop listening
    """

    def __init__(self, queue, *handlers):
        """Initialize the QueueListener.

        Args:
            queue: A multiprocessing-safe queue (e.g., from multiprocessing.Manager)
            *handlers: Variable number of logging handlers to forward records to
        """
        super().__init__()
        self.queue = queue
        self.handlers = handlers
        self._stop_event = threading.Event()

    def run(self):
        """Main loop that listens for log records and forwards them to handlers.

        Continuously polls the queue for log records and forwards them to
        appropriate handlers based on log level. Stops when a sentinel value
        (None) is received or when stop() is called.
        """
        while not self._stop_event.is_set():
            try:
                record = self.queue.get(timeout=1.0)
                if record is None:  # Sentinel to stop
                    break
                for handler in self.handlers:
                    # Only forward to handler if the record level meets the handler's level
                    if record.levelno >= handler.level:
                        handler.handle(record)
            except Exception:
                continue

    def stop(self):
        """Stop the listener thread.

        Sets the stop event and sends a sentinel value to the queue to ensure
        the thread exits cleanly.
        """
        self._stop_event.set()
        self.queue.put(None)  # Sentinel to stop the thread

    def flush(self, timeout=5.0):
        """Wait for all queued log records to be processed.
        
        Args:
            timeout (float): Maximum time to wait in seconds
            
        Returns:
            bool: True if queue was flushed successfully, False if timeout occurred
        """
        start_time = time.time()
        while not self.queue.empty() and (time.time() - start_time) < timeout:
            time.sleep(0.01)  # Small sleep to avoid busy waiting
        
        return self.queue.empty()


def setup_queue_logger(name, log_file, log_directory=None, console_level=logging.INFO,
                       file_level=logging.DEBUG, logger_object_level=logging.DEBUG):
    """Set up a queue-based logger for multiprocessing-safe logging.

    Creates a logger that uses a queue-based system for thread-safe logging
    in multiprocessing environments. The logger sends records to a shared
    queue, which are then consumed by a listener thread that forwards them
    to file and console handlers.

    Args:
        name (str): Name for the logger instance
        log_file (str): Filename for log file
        log_directory (str, optional): Directory for log files. If None, uses current directory
        console_level: Logging level for console output. Defaults to logging.INFO
        file_level: Logging level for file output. Defaults to logging.DEBUG
        logger_object_level: Logging level for logger object. Defaults to logging.DEBUG

    Returns:
        tuple: A tuple containing:
            - logger: Configured logger instance
            - queue: Multiprocessing-safe queue for log records
            - listener: Thread that consumes and forwards log records
    """
    # Get date and time string for log filename
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Create log directory if supplied
    if log_directory:
        if not os.path.exists(log_directory):
            os.makedirs(log_directory)
        log_file_name = f'{log_directory}/{log_file}_{date_and_time}.log'
    else:
        log_file_name = f'{log_file}_{date_and_time}.log'

    # Create handlers
    file_handler = logging.FileHandler(log_file_name, mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter(
        '%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_format)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Create manager and shared queue
    manager = Manager()
    queue = manager.Queue()
    queue_handler = QueueHandler(queue)

    # Set up logger
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logger_object_level)
    logger_object.addHandler(queue_handler)

    # Create and start listener
    listener = QueueListener(queue, file_handler, console_handler)
    listener.start()

    return logger_object, queue, listener


def setup_worker_logger(name, queue):
    """Set up a logger for worker processes that sends records to a shared queue.

    Creates a logger specifically designed for worker processes in multiprocessing
    environments. This logger sends all log records to a shared queue instead of
    directly to handlers, ensuring thread-safe logging across process boundaries.

    Args:
        name (str): Logger name (typically __name__ from the calling module)
        queue: Multiprocessing-safe queue for sending log records

    Returns:
        logging.Logger: Configured logger that sends records to the shared queue
    """
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logging.DEBUG)

    # Remove any existing handlers to avoid duplicates
    for handler in logger_object.handlers[:]:
        logger_object.removeHandler(handler)

    # Add queue handler
    queue_handler = QueueHandler(queue)
    logger_object.addHandler(queue_handler)

    return logger_object


class TqdmLogger:
    """File-like class redirecting tqdm progress bar to given logging logger."""
    def __init__(self, logger: logging.Logger):
        self.logger = logger

    def write(self, msg: str) -> None:
        # Remove carriage returns and empty lines from tqdm output
        if msg.strip() and msg.strip() != '\r':
            self.logger.info(msg.lstrip('\r'))

    def flush(self) -> None:
        pass


def createfolder(directory):
    """Creates a directory if it doesn't exist.

    Args:
        directory (str): Path to the directory to create.

    Returns:
        str: Path to the created directory.

    Raises:
        SystemExit: If the directory cannot be created due to an OSError.
    """

    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory
    except OSError:
        print(f'{"[ERROR]:":10} Error creating directory: {directory}')
        sys.exit(1)


def cprofile_to_csv(profile_binary_file):
    """Converts a cProfile.Profile object to CSV format.

    Takes a cProfile.Profile object, converts it to human-readable format, and 
    returns the data in CSV format for writing to file.

    Args:
        profile_binary_file (cProfile.Profile): A cProfile.Profile object.

    Returns:
        str: Human-readable data from the cProfile run in CSV format.

    Note:
        Adapted from: https://gist.github.com/ralfstx/a173a7e4c37afa105a66f371a09aa83e
    """
    out_stream = io.StringIO()
    pstats.Stats(profile_binary_file, stream=out_stream).sort_stats('cumtime').print_stats()
    result = out_stream.getvalue()
    
    # Remove header lines and keep only the data
    result = 'ncalls' + result.split('ncalls')[-1]
    
    # Convert each line to CSV format
    lines = [','.join(line.rstrip().split(None, 5)) for line in result.split('\n')]

    return '    \n'.join(lines)


def print_arguments(args, logger, __version__):
    """Prints the arguments to screen and log.

    Args:
        args (argparse.Namespace): Parsed command line arguments containing input file paths
            and other configuration options.
        logger (logging.Logger): Logger object.
        __version__ (str): Version of the script.
    """

    logger.info(f'{"[INFO]:":10} Plastid Annotation Validator (PAV) version {__version__} subcommand '
                f'`{args.subcommand_name}` was called with these arguments:\n')

    for parameter, value in args.__dict__.items():
        if parameter not in ['func', 'logger', 'log_directory', 'report_directory', 'subcommand_name']:
            logger.info(f'{" " * 10} {parameter}: {value}')
    logger.info('')


def check_dependencies(logger, entry='main'):
    """Checks for the presence of required external executables.

    Args:
        Nil.

    Returns:
        bool: True if all dependencies are found, False otherwise.
    """
    executables = ['mafft', 'trimal', 'blastn', 'julia']

    if entry == 'check':
        executables = ['mafft', 'trimal', 'blastn']
    elif entry == 'db_for_intgen':
        executables = ['makeblastdb']

    logger.info(f'{"[INFO]:":10} Checking for external dependencies:\n')

    all_executables_found = True
    for executable in executables:
        executable_loc = shutil.which(executable)
        if executable_loc:
            logger.info(f'{"":15} {executable:20} found at {executable_loc}')
        else:
            logger.info(f'{"":15} {executable:20} not found in your $PATH!')
            all_executables_found = False

    logger.info('')

    if all_executables_found:
        logger.info(f'{"":15} All external dependencies found!\n')
    else:
        logger.error(f'{"[ERROR]:":15} One or more dependencies not found!')
        exit_program()


def pad_seq(sequence):
    """Pads a sequence to a multiple of 3 with 'N' nucleotides.
    
    Takes a Bio.SeqRecord.SeqRecord object and ensures its sequence length is a 
    multiple of 3 by adding 'N' characters at the end if necessary. This is useful
    for ensuring DNA sequences can be properly translated.

    Args:
        sequence (Bio.SeqRecord.SeqRecord): Sequence record to pad.

    Returns:
        tuple: (Bio.SeqRecord.SeqRecord, bool) where the first element is the 
               padded sequence and the second indicates whether padding was applied.
    """
    remainder = len(sequence.seq) % 3
    
    if remainder == 0:
        return sequence, False
    else:
        # Add N's to make length a multiple of 3
        padding_length = 3 - remainder
        sequence.seq = sequence.seq + Seq('N' * padding_length)
        return sequence, True
    

def file_exists_and_not_empty(file_name):
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes
    """
    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def format_elapsed_time(elapsed_time):
    """
    Format elapsed time in seconds to hours, minutes, and seconds format.
    
    Args:
        elapsed_time (float): Elapsed time in seconds
        
    Returns:
        tuple: (hours, minutes, seconds, formatted_string)
    """
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)
    formatted_string = f'{hours}h {minutes}m {seconds}s ({elapsed_time:.2f} seconds)'
    
    return hours, minutes, seconds, formatted_string


def log_completion_time(start_time, logger=None, label='Completed'):
    """
    Log or print a standardized completion-time message given a start time.
    
    Args:
        start_time (float): Start timestamp from time.time().
        logger (logging.Logger, optional): Logger to write to; falls back to print if None.
        label (str): Message label preceding the duration (default: 'Completed').
    """
    try:
        elapsed_seconds = time.time() - start_time
        hours, minutes, seconds, _ = format_elapsed_time(elapsed_seconds)
        duration_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
        message = f"{"[INFO]:":10} {label} in {duration_str} ({elapsed_seconds:.2f} seconds)"
        if logger:
            logger.info(message)
        else:
            # If no logger, print a simplified message without the fixed-width tag
            print(f"{label} in {duration_str} ({elapsed_seconds:.2f} seconds)")
    except Exception:
        # Swallow any timing/logging errors to avoid masking real failures
        pass

def handle_process_result(result, id, operation_name, logger, log_listener):
    """Handle the result from a worker process, checking for errors and logging appropriately.
    
    Args:
        result: The result from the worker process (expected to be gene_id, (e, traceback.format_exc()) for errors)
        id (str): The id identifier for error reporting
        operation_name (str): The name of the operation (e.g., "MAFFT", "HMM generation")
        logger (logging.Logger): Logger object for output
        log_listener (QueueListener): QueueListener object for logging
        
    Returns:
        bool: True if result is successful, False if it's an error (script will exit)
        
    Raises:
        SystemExit: If any error is detected
    """

    # Check if result is an error tuple: (e, traceback.format_exc())
    if isinstance(result, tuple) and len(result) == 2 and isinstance(result[0], Exception):
        e, tb = result
        
        if isinstance(e, subprocess.CalledProcessError):
            logger.error(f'{"[ERROR]:":15} {operation_name} failed for ID {id}')
            logger.error(f'{"":15} Command: {" ".join(e.cmd)}')
            logger.error(f'{"":15} Return code: {e.returncode}')
            logger.error(f'{"":15} STDOUT: {e.stdout}')
            logger.error(f'{"":15} STDERR: {e.stderr}')
            logger.error(f'{"":15} Full traceback:')
            logger.error(tb)

            # Wait for log messages to be processed
            if log_listener:
                log_listener.flush(timeout=2.0)
                log_listener.stop()
                log_listener.join()
            sys.exit(1)
        else:
            logger.error(f'{"[ERROR]:":15} {operation_name} failed for ID {id}')
            logger.error(f'{"":15} Exception: {e}')
            logger.error(f'{"":15} Full traceback:')
            logger.error(tb)

            # Wait for log messages to be processed
            if log_listener:
                log_listener.flush(timeout=2.0)
                log_listener.stop()
                log_listener.join()
            sys.exit(1)
    else:
        logger.error(f'{"[ERROR]:":15} {operation_name} failed for ID {id} (unknown error)')

        # Wait for log messages to be processed
        if log_listener:
            log_listener.flush(timeout=2.0)
            log_listener.stop()
            log_listener.join()
        sys.exit(1)


def build_adjusted_genbank_io(gbk_file_path, fasta_file_path, logger=None):
    """Create an in-memory adjusted GenBank file object after fixing the LOCUS line.

    Reads the corresponding FASTA file to compute the total number of base pairs, then adjusts
    the LOCUS line in the GenBank file to reflect this length and enforces four spaces between
    the last three fields while preserving the original prefix spacing. Returns both a StringIO
    handle of the adjusted content and the adjusted content string for optional previewing.

    Args:
        gbk_file_path (str): Path to the GenBank file to adjust.
        fasta_file_path (str): Path to the FASTA file used to determine sequence length.
        logger (logging.Logger, optional): Logger instance for messages and warnings.

    Returns:
        tuple[io.StringIO, str]: (in_memory_gbk_io, adjusted_gbk_content)

    Raises:
        ValueError: If the GenBank file is empty.
    """
    # Read FASTA to get total bp length
    fasta_bp = 0
    try:
        for rec in SeqIO.parse(fasta_file_path, "fasta"):
            fasta_bp += len(rec.seq)
    except Exception:
        fasta_bp = 0

    # Read and adjust GenBank LOCUS line (preserve original spacing except for required tweaks)
    # Handle both gzipped and uncompressed files
    if gbk_file_path.endswith('.gz'):
        with gzip.open(gbk_file_path, 'rt') as gbk_handle:
            gbk_lines = gbk_handle.readlines()
    else:
        with open(gbk_file_path, 'r') as gbk_handle:
            gbk_lines = gbk_handle.readlines()
    if not gbk_lines:
        if logger:
            logger.warning(f"{'[WARNING]:':10} Empty GenBank file: {gbk_file_path}")
        raise ValueError(f"Empty GenBank file: {gbk_file_path}")

    original_first_line = gbk_lines[0].rstrip('\n')

    adjusted_first = original_first_line
    if original_first_line.startswith('LOCUS'):
        # Replace length just before 'bp' with FASTA-derived bp while keeping other spacing
        if fasta_bp > 0:
            adjusted_first = re.sub(r'0 bp', f'{fasta_bp} bp', adjusted_first)

        # Enforce four spaces between the last three fields only, preserving the prefix
        parts = adjusted_first.rsplit(None, 3)
        if len(parts) == 4:
            head, f3, f2, f1 = parts  # head contains original spacing for the prefix
            adjusted_first = f"{head}    {f3}    {f2}    {f1}"

    adjusted_first_line = adjusted_first + '\n'
    adjusted_gbk_content = adjusted_first_line + ''.join(gbk_lines[1:])
    gbk_io = io.StringIO(adjusted_gbk_content)

    return gbk_io, adjusted_gbk_content


def pad_seq(sequence):
    """Pads a sequence to a multiple of 3 with 'N' nucleotides.
    
    Takes a Bio.SeqRecord.SeqRecord object and ensures its sequence length is a 
    multiple of 3 by adding 'N' characters at the end if necessary. This is useful
    for ensuring DNA sequences can be properly translated.

    Args:
        sequence (Bio.SeqRecord.SeqRecord): Sequence record to pad.

    Returns:
        tuple: (Bio.SeqRecord.SeqRecord, bool) where the first element is the 
               padded sequence and the second indicates whether padding was applied.
    """
    remainder = len(sequence.seq) % 3
    
    if remainder == 0:
        return sequence, False
    else:
        # Add N's to make length a multiple of 3
        padding_length = 3 - remainder
        sequence.seq = sequence.seq + Seq('N' * padding_length)
        return sequence, True
    

def log_separator(logger, length=100):
    """
    Log a separator line with consistent formatting.
    
    Args:
        logger: Logger instance for logging messages
        length (int): Length of the separator line (default: 100)
    """
    logger.info(f"{"":10} {'-'*length}")


def parse_genbank_file(gbk_file_path, logger=None):
    """
    Open and parse a GenBank file (zipped or unzipped) with warning filtering.
    
    This function handles both compressed (.gz) and uncompressed GenBank files,
    automatically detecting the file type and applying appropriate warning filters
    to suppress Biopython parser warnings about malformed locus lines and
    sequence length mismatches.
    
    Args:
        gbk_file_path (str): Path to the GenBank file (.gb, .gbk, .gb.gz, .gbk.gz)
        logger (logging.Logger, optional): Logger instance for messages and warnings
        
    Returns:
        list: List of Bio.SeqRecord.SeqRecord objects from the GenBank file
        
    Raises:
        FileNotFoundError: If the GenBank file doesn't exist
        ValueError: If the GenBank file is empty or contains no valid records
        Exception: For other parsing errors
    """
    import warnings
    
    # Check if file exists
    if not os.path.exists(gbk_file_path):
        raise FileNotFoundError(f"GenBank file not found: {gbk_file_path}")
    
    try:
        # Parse GenBank file with warning filtering
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            warnings.filterwarnings("ignore", message=".*malformed locus line.*")
            warnings.filterwarnings("ignore", message=".*Expected sequence length.*")
            warnings.filterwarnings("ignore", message=".*Attempting to parse malformed locus line.*")
            warnings.filterwarnings("ignore", message=".*Some fields may be wrong.*")
            
            if gbk_file_path.endswith('.gz'):
                with gzip.open(gbk_file_path, 'rt') as handle:
                    records = list(SeqIO.parse(handle, 'genbank'))
            else:
                with open(gbk_file_path, 'r') as handle:
                    records = list(SeqIO.parse(handle, 'genbank'))
        
        if not records:
            if logger:
                logger.warning(f"{'[WARNING]:':10} No sequences found in GenBank file: {gbk_file_path}")
            raise ValueError(f"No sequences found in GenBank file: {gbk_file_path}")
        
        if logger:
            logger.debug(f"{'[DEBUG]:':10} Successfully parsed {len(records)} records from {gbk_file_path}")
        
        return records
        
    except Exception as e:
        if logger:
            logger.error(f"{'[ERROR]:':10} Error parsing GenBank file {gbk_file_path}: {e}")
        raise


def exit_program():
    """Exit the program with a clean shutdown.
    
    Args:
        None
    """
    log_manager.cleanup()
    sys.exit(1) 