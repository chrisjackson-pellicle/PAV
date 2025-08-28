#!/usr/bin/env python
"""
Plastid Annotation Validator (PAV) - Annotation and Validation Module

This module provides comprehensive functionality for annotating plastid genomes and validating
their gene annotations against reference data. It serves as the core processing engine for the
PAV pipeline, handling genome annotation, gene validation, alignment generation, and intergenic
region analysis.

Main Functions:
    - annotate_genomes(): Annotates plastid genomes using the chloë annotation pipeline
    - check_genes(): Validates gene annotations against reference median lengths and synonyms
    - align_genes(): Generates multiple sequence alignments for annotated genes
    - query_intergenic_regions(): Performs BLAST analysis of intergenic regions
    - convert_gbk_to_embl(): Converts GenBank files to EMBL format with locus tags

Key Features:
    - Multi-process parallel processing for improved performance
    - Comprehensive logging with queue-based multiprocessing support
    - Gene synonym mapping for standardized gene naming
    - Length validation against reference median values
    - Multiple alignment formats (CDS, rRNA, per-gene)
    - Intergenic region extraction and BLAST analysis
    - Robust error handling and cleanup

Dependencies:
    - Biopython: Sequence parsing and manipulation
    - chloë: Genome annotation pipeline
    - MAFFT: Multiple sequence alignment
    - BLAST+: Sequence similarity search
    - TrimAl: Alignment trimming

Usage:
    This module is typically called through the PAV command-line interface:
    $ pav annotate_and_check [options]

Author: Chris Jackson chris.jackson@rbg.vic.gov.au
License: See LICENSE file
"""

import os
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import subprocess
import tempfile
import glob
import pandas as pd
import importlib.resources
import time
import re

# Import PAV modules:
from plastid_annotation_validator.version import __version__
from plastid_annotation_validator import utils

# Initialise logger objects to None
logger = None
log_queue = None
log_listener = None

def load_gene_synonyms(data_dir_base=None):
    """
    Load gene synonyms from the package text file.
    
    Args:
        None

    Returns:
        dict: Dictionary mapping gene names to their standardized synonyms
        
    """

    # Require provided data_dir_base and do not fall back
    if data_dir_base is None:
        raise ValueError("data_dir_base must be provided to load_gene_synonyms()")
    candidate_path = os.path.join(data_dir_base, 'gene_synonyms.txt')
    if not os.path.exists(candidate_path):
        raise FileNotFoundError(f"Gene synonyms file not found at: {candidate_path}")
    
    gene_synonyms = {}
    with open(candidate_path, 'r') as synonyms_file:
        for line in synonyms_file:
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            # Parse the mapping: "GenBank_name -> Standard_name"
            if ' = ' in line:
                genbank_name, standard_name = line.split(' = ', 1)
                gene_synonyms[genbank_name.strip()] = standard_name.strip()

    if logger:
        logger.debug(f"{"[DEBUG]:":10} Successfully loaded {len(gene_synonyms)} gene synonyms from package resources")
    return gene_synonyms
        

def load_gene_median_lengths(data_dir_base=None):
    """ 
    Load gene median lengths from the package CSV file.
    
    Args:
        logger: Logger instance for logging messages

    Returns:
        dict: Dictionary mapping gene names to their median lengths
        
    """

    # Require provided data_dir_base and do not fall back
    if data_dir_base is None:
        raise ValueError("data_dir_base must be provided to load_gene_median_lengths()")

    candidate_path = os.path.join(data_dir_base, 'plDNA_genes_median_lengths.csv')
    if not os.path.exists(candidate_path):
        raise FileNotFoundError(f"Median lengths CSV not found at: {candidate_path}")

    with open(candidate_path, 'r') as csv_file:
        gene_lengths_df = pd.read_csv(csv_file)

    if logger:
        logger.debug(f"{"[DEBUG]:":10} Successfully loaded reference gene median lengths from package resources")
        logger.debug(f"{"[DEBUG]:":10} Loaded {len(gene_lengths_df)} gene entries")
    
    # Convert to dictionary for easier access
    gene_median_lengths = dict(zip(gene_lengths_df['Key'], gene_lengths_df['Value']))
    if logger:
        logger.debug(f"{"[DEBUG]:":10} Loaded {len(gene_median_lengths)} gene median lengths")
    
    
    return gene_median_lengths


def parse_gbk_genes(gbk_file_path, logger, gene_synonyms=None):
    """
    Parse gene annotations from a GenBank file and extract gene lengths, CDS, rRNA, and tRNA information.
    
    This function processes a GenBank file to extract comprehensive gene annotation information:
    1. Parses CDS, rRNA, and tRNA features from the GenBank file
    2. Extracts gene names and applies synonym mapping for standardization
    3. Creates individual entries for each gene copy with appropriate suffixes
    4. Extracts sequence information for translation validation
    5. Organizes information into structured dictionaries for downstream analysis
    
    Args:
        gbk_file_path (str): Path to the GenBank file (.gb or .gbk)
        logger: Logger instance for logging messages and debugging
        gene_synonyms (dict, optional): Dictionary mapping gene names to standardized synonyms.
            Used to normalize gene naming across different annotation sources.
        
    Returns:
        tuple: (gene_lengths, gene_cds_info, gene_rRNA_info, gene_tRNA_info) where:
            - gene_lengths (dict): Dictionary mapping gene names to length information containing:
                - 'annotated_gene_length': Individual gene length for this copy
                - 'copies': Always 1 (each copy is stored separately)
                - 'copy_lengths': List containing single length for this copy
            - gene_cds_info (dict): Dictionary mapping gene names to CDS information for translation checking
            - gene_rRNA_info (dict): Dictionary mapping gene names to rRNA information
            - gene_tRNA_info (dict): Dictionary mapping gene names to tRNA information
            
    Note:
        Multiple copies of the same gene are stored as separate entries with suffixes:
        - First copy: original gene name (e.g., 'rps12')
        - Subsequent copies: gene name with '_copy_N' suffix (e.g., 'rps12_copy_2', 'rps12_copy_3')
    
    Raises:
        ValueError: If no records are found in the GenBank file or if a feature lacks a gene name
        Exception: Various exceptions that may occur during file parsing or sequence extraction
        
    Note:
        The function handles multiple copies of the same gene by creating separate entries
        for each copy with appropriate suffixes. Gene names are standardized using the
        provided synonym mapping if available.
    """
    gene_lengths = {}
    gene_cds_info = {}
    gene_rRNA_info = {}
    gene_tRNA_info = {}
    
    try:
        records = utils.parse_genbank_file(gbk_file_path, logger)
        
        # Process each record (should be only one as I've split multi-fasta in to individual files)
        for record in records:
            # Initialize data structures for this record
            gene_lengths = defaultdict(list)
            gene_data = dict()
            gene_cds_locations = defaultdict(list)
            gene_rRNA_locations = defaultdict(list)
            gene_tRNA_locations = defaultdict(list)
            
            # Process each feature in the record
            for feature in record.features:
                if feature.type not in ["CDS", "rRNA", "tRNA"]:
                    continue
                    
                # Extract gene name from qualifiers
                gene_name = extract_gene_name(feature)
                if not gene_name:
                    raise ValueError(f"No gene name found for feature: {feature}")
                
                # Apply gene synonym mapping if available
                mapped_gene_name = gene_name
                if gene_synonyms and gene_name in gene_synonyms:
                    mapped_gene_name = gene_synonyms[gene_name]
                    logger.debug(f"{"[DEBUG]:":10} Mapped gene name: {gene_name} -> {mapped_gene_name}")

                # Calculate gene length and extract sequence
                gene_length = len(feature.location)
                extracted_seq = feature.location.extract(record.seq)
                
                # Check for length mismatch between feature location and extracted sequence
                length_mismatch_msg = None
                if gene_length != len(extracted_seq):
                    logger.debug(f"{"[DEBUG]:":10} Gene {mapped_gene_name} length mismatch: Expected {gene_length} "
                                 f"!= Extracted {len(extracted_seq)}")
                    length_mismatch_msg = (f"Gene {mapped_gene_name} length mismatch: Expected {gene_length} "
                                           f"!= Extracted {len(extracted_seq)}")
                
                # Store gene length information
                gene_lengths[mapped_gene_name].append(gene_length)
                
                # Store feature-specific information
                feature_info = {
                    'location': feature.location,
                    'strand': feature.location.strand,
                    'length': gene_length,
                    'length_mismatch_msg': length_mismatch_msg
                }
                
                if feature.type == "CDS":
                    feature_info['cds_seq'] = extracted_seq
                    gene_cds_locations[mapped_gene_name].append(feature_info)
                elif feature.type == "rRNA":
                    feature_info['rRNA_seq'] = extracted_seq
                    gene_rRNA_locations[mapped_gene_name].append(feature_info)
                elif feature.type == "tRNA":
                    feature_info['tRNA_seq'] = extracted_seq
                    gene_tRNA_locations[mapped_gene_name].append(feature_info)
            
            # Process each gene and create individual entries for each copy
            for gene_name, lengths in gene_lengths.items():
                # Get the feature lists for this gene
                cds_features = gene_cds_locations.get(gene_name, [])
                rrna_features = gene_rRNA_locations.get(gene_name, [])
                trna_features = gene_tRNA_locations.get(gene_name, [])
                
                # Determine how many copies we have (use the maximum from any feature type)
                total_copies = max(len(lengths), len(cds_features), len(rrna_features), len(trna_features))
                
                if total_copies == 1:
                    # Single copy - use original gene name
                    copy_name = gene_name
                    
                    # Store gene length information
                    gene_data[copy_name] = {
                        'annotated_gene_length': lengths[0],
                        'copies': 1,
                        'copy_lengths': lengths
                    }
                    
                    # Store feature information
                    if cds_features:
                        gene_cds_info[copy_name] = cds_features
                    if rrna_features:
                        gene_rRNA_info[copy_name] = rrna_features
                    if trna_features:
                        gene_tRNA_info[copy_name] = trna_features
                        
                else:
                    # Multiple copies - create individual entries with suffixes
                    logger.debug(f"{"[DEBUG]:":10} Gene {gene_name} has {total_copies} copies, creating individual entries")
                    
                    for copy_num in range(1, total_copies + 1):
                        copy_name = f"{gene_name}_copy_{copy_num}"
                        
                        # Get the length for this copy 
                        copy_length = lengths[copy_num - 1]
                        
                        # Store gene length information
                        gene_data[copy_name] = {
                            'annotated_gene_length': copy_length,
                            'copies': len(lengths),
                            'copy_lengths': lengths
                        }
                        
                        # Store feature information for this copy
                        if copy_num <= len(cds_features):
                            gene_cds_info[copy_name] = [cds_features[copy_num - 1]]
                        if copy_num <= len(rrna_features):
                            gene_rRNA_info[copy_name] = [rrna_features[copy_num - 1]]
                        if copy_num <= len(trna_features):
                            gene_tRNA_info[copy_name] = [trna_features[copy_num - 1]]
        
        logger.debug(f"{"[DEBUG]:":10} Parsed {len(gene_lengths)} genes (CDS, rRNA, and tRNA) from GenBank file")
        return gene_data, gene_cds_info, gene_rRNA_info, gene_tRNA_info
        
    except Exception as e:
        logger.error(f"{"[ERROR]:":10} Error parsing GenBank file: {e}")
        logger.error(traceback.format_exc())
        raise


def extract_gene_name(feature):
    """
    Extract gene name from a Biopython feature's qualifiers.
    
    Args:
        feature: Biopython SeqFeature object
        
    Returns:
        str: Gene name if found, None otherwise
    """
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    return None


def check_gene_translation(gene_name, cds_info_list, logger):
    """
    Check if a gene's CDS can be translated correctly.
    
    Args:
        gene_name (str): Name of the gene
        cds_info_list (list): List of CDS information dictionaries containing 'cds_seqrecord'
        logger: Logger instance for logging messages
        
    Returns:
        list: List of warning dictionaries for translation issues
    """
    warnings = []

    
    try:
        for i, cds_info in enumerate(cds_info_list):  # can be multiple CDSs per gene due to multiple copies

            # Check for length mismatch
            if cds_info['length_mismatch_msg']:
                warnings.append(f"{cds_info['length_mismatch_msg']} for copy {i+1} - check this gene manually!")

            cds_seq = cds_info['cds_seq']  # Already extracted sequence

            # Create SeqRecord object
            cds_seqrecord = SeqRecord(cds_seq, id=f"{gene_name}_copy_{i+1}", description="")
            
            # Check if sequence length is a multiple of 3
            if len(cds_seqrecord.seq) % 3 != 0:
                warnings.append(f"CDS length {len(cds_seqrecord.seq)} is not a multiple of 3 for copy {i+1}")
                
            # Check for start codon (ATG) at the beginning
            start_codon = str(cds_seqrecord.seq[:3])

            if start_codon != "ATG":
                # warnings.append(f"Start codon is {start_codon}, expected ATG for copy {i+1}")
                warnings.append(f"Start codon is {start_codon}, expected ATG")
            
            # Check for stop codon at the end
            stop_codon = str(cds_seqrecord.seq[-3:])

            stop_codons = ["TAA", "TAG", "TGA"]
            if stop_codon not in stop_codons:
                warnings.append(f"Stop codon is {stop_codon}, expected one of {stop_codons} for copy {i+1}")

            # Try to translate the sequence
            try:
                # Pad sequence to multiple of 3
                cds_seqrecord, _ = utils.pad_seq(cds_seqrecord)
                translated = cds_seqrecord.translate()

                # Check for premature stop codons (* in translation)
                if "*" in str(translated.seq)[:-1]:  # Exclude the last position which should be a stop codon
                    warnings.append(f"Premature stop codon(s) found in translation for copy {i+1}")
                
            except Exception as e:
                warnings.append(f"Translation failed: {e} for copy {i+1}")
        
        return warnings
        
    except Exception as e:
        logger.error(f"{"[ERROR]:":10} Error checking translation for gene {gene_name}: {e}")
        raise e


def check_single_sample_genes(sample_data, gene_median_lengths, min_threshold, max_threshold, gene_synonyms=None,
                              log_queue=None):
    """
    Validate gene annotations for a single plastid genome sample by checking gene lengths and translation quality.
    
    This function performs comprehensive validation of gene annotations for one sample by:
    1. Parsing the GenBank files to extract gene annotations and sequences (supports multiple files per sample)
    2. Comparing gene lengths against reference median lengths from a curated database
    3. Checking translation quality for CDS genes to identify potential annotation issues
    4. Applying gene synonym mapping for standardized gene naming
    5. Generating detailed validation results for each gene
    
    Args:
        sample_data (dict): Dictionary containing sample information with keys:
            - 'sample_name': Name of the sample being processed
            - 'gbk_files': List of paths to annotated GenBank files
        gene_median_lengths (dict): Dictionary mapping gene names to their reference median lengths
        min_threshold (float): Minimum acceptable length as percentage of median (e.g., 0.8 for 80%)
        max_threshold (float): Maximum acceptable length as percentage of median (e.g., 1.2 for 120%)
        gene_synonyms (dict, optional): Dictionary mapping gene names to standardized synonyms.
            Used to normalize gene naming across different annotation sources.
        log_queue (queue.Queue, optional): Multiprocessing-safe queue for logging messages.
            Required for proper logging in multiprocessing environments.
    
    Returns:
        tuple: (success, result) where:
            - success (bool): True if processing completed successfully, False if an error occurred
            - result: If success=True, a dictionary containing:
                - 'sample_name': Name of the processed sample
                - 'sequence_gene_data': Dictionary mapping sequence names to gene data containing:
                    - 'gbk_file': Path to GenBank file for this sequence
                    - 'gene_cds_info': CDS information for this sequence
                    - 'gene_rRNA_info': rRNA information for this sequence
                    - 'gene_tRNA_info': tRNA information for this sequence
                    - 'sequence_gene_data': Dictionary mapping gene names to validation results containing:
                        - 'gene_name': Gene identifier
                        - 'annotated_gene_length': Observed gene length in base pairs
                        - 'copies': Number of copies of this gene in the genome
                        - 'copy_lengths': List of lengths for each copy
                        - 'threshold_min': Minimum threshold used for validation
                        - 'threshold_max': Maximum threshold used for validation
                        - 'median_length': Reference median length (if available)
                        - 'issue': Validation status ('OK', 'too_short', 'too_long', 'not_in_median_lengths')
                        - 'details': Detailed description of any issues found
                        - 'translation_status': Translation validation status ('OK' or 'FAIL')
                        - 'translation_details': Details about translation issues (if any)
                - 'total_genes': Total number of genes processed
                - 'genes_with_warnings': Count of genes that failed validation
                - 'gene_cds_info': Combined CDS location information for translation checking
                - 'gene_rRNA_info': Combined rRNA location information
                - 'gene_tRNA_info': Combined tRNA location information
            If success=False, result is a tuple (exception, traceback) for error handling
    
    Raises:
        ValueError: If there are inconsistencies between gene lengths and CDS information
        Exception: Various exceptions that may occur during file parsing or processing
        
    Note:
        This function is designed to be called by multiprocessing workers and includes
        comprehensive error handling to ensure robust processing of individual samples.
        It sets up a worker-specific logger to ensure proper logging in multiprocessing
        environments.
    """
    try:
        worker_logger = utils.setup_worker_logger(__name__, log_queue)

        sample_name = sample_data.get('sample_name', 'unknown')
        gbk_files = sample_data['gbk_files']

        length_warnings = set()
        warnings_to_log_to_screen = set()
        
        # Process each file for a sample and combine the results
        sequence_gene_data = {}
        combined_gene_cds_info = {}
        combined_gene_rRNA_info = {}
        combined_gene_tRNA_info = {}
        
        # Process each file for a sample and combine the results
        for gbk_file in gbk_files:
            # Parse gene lengths and CDS information from GenBank file (adjusted in-memory)
            gene_data, gene_cds_info, gene_rRNA_info, gene_tRNA_info = parse_gbk_genes(gbk_file,
                                                                                       worker_logger, gene_synonyms)

            # Store individual sequence data for reporting
            seq_filename = os.path.basename(gbk_file)
            # Extract sequence name by removing file extensions
            seq_name = seq_filename
            for ext in ['.chloe.gbk', '.round1.chloe.gbk', 'round2.chloe.gbk']:
                if seq_name.endswith(ext):
                    seq_name = seq_name[:-len(ext)]
                    break

            # Store individual sequence data with validation fields initialized
            sequence_gene_data[seq_name] = {
                'gbk_file': gbk_file,
                'gene_cds_info': gene_cds_info,
                'gene_rRNA_info': gene_rRNA_info,
                'gene_tRNA_info': gene_tRNA_info,
                'sequence_gene_data': {}
            }
            
            # Process genes for this sequence and add validation fields
            for gene_name, gene_info in gene_data.items():
                
                # Add validation fields to gene info
                gene_info_with_validation = {
                    'threshold_min': min_threshold,
                    'threshold_max': max_threshold,
                    'ref_median_length': None,
                    'issue': None,
                    'details': None,
                    'translation_status': None,
                    'translation_details': None,
                }
                gene_info_with_validation.update(gene_info)

                sequence_gene_data[seq_name]['sequence_gene_data'][gene_name] = gene_info_with_validation

            # Combine CDS, rRNA, and tRNA info from each sequence into a single dict
            for gene_name, info in gene_cds_info.items():
                if gene_name in combined_gene_cds_info:
                    combined_gene_cds_info[gene_name].extend(info)
                else:
                    combined_gene_cds_info[gene_name] = info.copy()
            
            for gene_name, info in gene_rRNA_info.items():
                if gene_name in combined_gene_rRNA_info:
                    combined_gene_rRNA_info[gene_name].extend(info)
                else:
                    combined_gene_rRNA_info[gene_name] = info.copy()
            
            for gene_name, info in gene_tRNA_info.items():
                if gene_name in combined_gene_tRNA_info:
                    combined_gene_tRNA_info[gene_name].extend(info)
                else:
                    combined_gene_tRNA_info[gene_name] = info.copy()
        
        # Get all unique gene names from all sequences
        all_gene_names = set()
        for seq_data in sequence_gene_data.values():
            all_gene_names.update(seq_data['sequence_gene_data'].keys())

        all_gene_names_no_copy_suffix = set([gene_name.split('_copy_')[0] for gene_name in all_gene_names])
        
        # Check for genes in CDS info that are not in gene lengths
        cds_only_genes = set(combined_gene_cds_info.keys()) - all_gene_names
        if cds_only_genes:
            return False, (ValueError(f"Genes in CDS info but not in gene lengths: {list(cds_only_genes)[:10]}"), traceback.format_exc())
        
        # Check for genes present in gene_lengths but not in gene_median_lengths, and vice versa
        genes_in_lengths_not_in_median = list(all_gene_names_no_copy_suffix - set(gene_median_lengths.keys()))
        genes_in_median_not_in_lengths = list(set(gene_median_lengths.keys()) - all_gene_names_no_copy_suffix)
        
        if genes_in_lengths_not_in_median:
            worker_logger.debug(f"{"[DEBUG]:":10} Genes in sample {sample_name} but not in median lengths: "
                                f"{genes_in_lengths_not_in_median}...")
        if genes_in_median_not_in_lengths:
            worker_logger.debug(f"{"[DEBUG]:":10} Genes in median lengths but not in sample {sample_name}: "
                                f"{genes_in_median_not_in_lengths}...")
       
        # Check each gene copy against median lengths and translation
        for seq_name, seq_data in sequence_gene_data.items():
            gene_info_dict = seq_data['sequence_gene_data']

            for gene_name, gene_info in gene_info_dict.items():
                annotated_gene_length = gene_info['annotated_gene_length']
                gene_name_no_copy_suffix = gene_name.split('_copy_')[0]
                
                # Check length against median
                if gene_name_no_copy_suffix in gene_median_lengths:
                    median_length = gene_median_lengths[gene_name_no_copy_suffix]
                    min_expected = median_length * min_threshold
                    max_expected = median_length * max_threshold

                    # Add median length info to gene_info
                    gene_info['ref_median_length'] = median_length
       
                    # Check if gene is too short
                    if annotated_gene_length < min_expected:
                        warning_msg = (f"Too short: {annotated_gene_length} bp < {min_expected:.0f} bp "
                                       f"({min_threshold*100}% of ref median)")
                        gene_info['issue'] = 'too_short'
                        gene_info['details'] = warning_msg
                        length_warnings.add(gene_name)

                    # Check if gene is too long
                    elif annotated_gene_length > max_expected:
                        warning_msg = (f"Too long: {annotated_gene_length} bp > {max_expected:.0f} bp "
                                       f"({max_threshold*100}% of median)")
                        gene_info['issue'] = 'too_long'
                        gene_info['details'] = warning_msg
                        length_warnings.add(gene_name)

                    else:
                        gene_info['issue'] = 'OK'
                        gene_info['details'] = 'OK'

                else:
                    worker_logger.debug(f"{"[DEBUG]:":10} Gene {gene_name_no_copy_suffix} from sample {sample_name} not in gene_median_lengths")
                    warnings_to_log_to_screen.add(f"Gene {gene_name_no_copy_suffix} from sample {sample_name} not in gene_median_lengths")
                    gene_info['issue'] = 'not_in_median_lengths'
                    gene_info['details'] = 'Not in reference median lengths'
                    length_warnings.add(gene_name)
                
                # Check translation for CDS genes
                if gene_name in combined_gene_cds_info:
                    translation_issues = check_gene_translation(gene_name, combined_gene_cds_info[gene_name], logger)

                    if translation_issues:
                        gene_info['translation_status'] = 'FAIL'
                        gene_info['translation_details'] = translation_issues
                        length_warnings.add(gene_name)
                    else:
                        gene_info['translation_status'] = 'OK'
                        gene_info['translation_details'] = 'OK'

        # Collate results for return
        result_dict = {
            'sample_name': sample_name,
            'total_genes': len(all_gene_names),
            'genes_with_warnings': len(length_warnings),
            'gene_cds_info': combined_gene_cds_info,
            'gene_rRNA_info': combined_gene_rRNA_info,
            'gene_tRNA_info': combined_gene_tRNA_info,
            'missing_genes': genes_in_median_not_in_lengths,
            'sequence_gene_data': sequence_gene_data,
            'warnings_to_log_to_screen': warnings_to_log_to_screen
        }

        return True, result_dict
        
    except Exception as e:
        return False, (e, traceback.format_exc())


def check_genes(gene_median_lengths, annotated_genomes_dict, min_threshold, max_threshold, 
                report_directory, log_queue, pool_size=1,  gene_synonyms=None):
    """
    Validate gene annotations across multiple plastid genomes by comparing gene lengths against reference median values.
    
    This function processes annotated plastid genomes and validates gene annotations by:
    1. Parsing GenBank files to extract gene annotations and sequences (supports multiple files per sample)
    2. Comparing gene lengths against reference median lengths from a curated database
    3. Applying gene synonym mapping for standardized gene naming
    4. Generating comprehensive reports of validation results
    5. Using multiprocessing for efficient parallel processing
    
    Args:
        gene_median_lengths (dict): Dictionary mapping gene names to their reference median lengths
        annotated_genomes_dict (dict): Dictionary mapping sample names to dicts containing:
            - 'gbk': Path to annotated GenBank file (single) or list of paths (multi-sequence)
            - 'is_multi_sequence': Boolean indicating if this is a multi-sequence sample
            - 'sequence_count': Number of sequences in this sample
        min_threshold (float): Minimum acceptable length as percentage of median (e.g., 80.0 for 80%)
        max_threshold (float): Maximum acceptable length as percentage of median (e.g., 120.0 for 120%)
        report_directory (str): Directory path where validation reports will be written
        log_queue (queue.Queue): Multiprocessing-safe queue for logging messages
        pool_size (int, optional): Number of worker processes for parallel processing. Defaults to 1.
        gene_synonyms (dict, optional): Dictionary mapping gene names to standardized synonyms.
            Used to normalize gene naming across different annotation sources.

    Returns:
        dict: Dictionary mapping sample names to validation results containing:
            - 'sequence_gene_data': Dictionary mapping sequence names to gene data
            - 'gene_cds_info': Combined CDS information across all sequences
            - 'gene_rRNA_info': Combined rRNA information across all sequences
            - 'gene_tRNA_info': Combined tRNA information across all sequences
            - 'genes_with_warnings': Count of genes that failed validation
            - 'total_genes': Total number of genes processed
            - 'error': Error message if processing failed (optional)
    
    Raises:
        SystemExit: If any sample fails to process and error handling is triggered
        
    Note:
        The function uses multiprocessing to process samples in parallel, with each sample
        being validated independently. For multi-sequence samples, all sequences are processed
        together and gene data is collated across all sequences. Results are collected and 
        combined into comprehensive reports including both individual sample reports and a 
        combined summary.
    """

    all_sample_results = {}
    total_warnings = 0
    
    # Prepare data for multiprocessing
    sample_data_list = []
    for sample_name, data in annotated_genomes_dict.items():
        # Handle both single files and lists of files (for multi-sequence samples)
        if isinstance(data['gbk'], list):
            # Multi-sequence sample - collect all files
            gbk_files = data['gbk']
        else:
            # Single sequence sample - convert to lists for consistent processing
            gbk_files = [data['gbk']]
        
        sample_data = {
            'gbk_files': gbk_files,
            'sample_name': sample_name
        }
        sample_data_list.append(sample_data)

    # Process samples with multiprocessing and progress bar

    warnings_to_log_to_screen = []
    with ProcessPoolExecutor(max_workers=pool_size) as executor:
        # Submit all tasks
        future_to_sample = {
            executor.submit(check_single_sample_genes, 
                            sample_data, 
                            gene_median_lengths,
                            min_threshold, 
                            max_threshold, 
                            gene_synonyms,
                            log_queue):  sample_data['sample_name'] 
                            
                            for sample_data in sample_data_list
        }
        
        # Process results with progress bar
        for future in tqdm(as_completed(future_to_sample), total=len(future_to_sample),
                           desc=f"{"[INFO]:":10} {"Validating genes":<20}", file=sys.stdout):
            sample_name = future_to_sample[future]

            success, result = future.result()
            
            if success:

                all_sample_results[sample_name] = result
                total_warnings += result['genes_with_warnings']
                warnings_to_log_to_screen.extend(result['warnings_to_log_to_screen'])

            else:
                utils.log_manager.handle_error(result[0], result[1], "check_genes()", sample_name)

    for warning in warnings_to_log_to_screen:
        logger.warning(f"{"[WARNING]:":10} {warning}")

    # Generate and write report
    write_gene_length_report(all_sample_results, logger, min_threshold, max_threshold, report_directory,
                             gene_median_lengths)
    
    utils.log_separator(logger)

    return all_sample_results


def write_gene_length_report(all_results, logger, min_threshold, max_threshold, report_directory, gene_median_lengths):
    """
    Generate comprehensive TSV reports for gene validation results across all samples.
    
    This function creates detailed reports of gene validation results including:
    1. Individual sample reports with gene-by-gene validation details
    2. Combined report with all samples and genes in a single file
    3. Summary statistics for warnings and validation issues
    4. Detailed information about gene lengths, copies, and translation quality
    
    Args:
        all_results (dict): Dictionary mapping sample names to validation results
        logger: Logger instance for logging messages and progress updates
        min_threshold (float): Minimum acceptable length as percentage of median
        max_threshold (float): Maximum acceptable length as percentage of median
        report_directory (str): Directory path where all reports will be written
        gene_median_lengths (dict): Dictionary mapping gene names to their reference median lengths
    
    Returns:
        None: Reports are written directly to files in the specified directory
    """
    
    def create_tsv_header():
        """Create TSV header with sequence column."""
        
        return ["Sequence", "Gene_Name", "Annotation_Length", "Median_Length_Refs", 
                "Min_Threshold", "Max_Threshold", "Copies", 
                "Warning", "Translation_Details"]
    
    def format_tsv_row(sequence_name, gene_name, annotation_length, median_length, 
                      min_threshold, max_threshold, copies, 
                      warning_msg, translation_details):
        """Format a single TSV row with sequence column."""

        median_str = str(median_length) if median_length is not None else 'NA'
        
        return [str(sequence_name), str(gene_name), str(annotation_length), median_str, 
                str(min_threshold), str(max_threshold), str(copies), 
                str(warning_msg), str(translation_details)]
    
    def process_gene_data(gene_name, gene_info, sequence_name):
        """Process gene data and return formatted row data."""
        # Extract data from gene_info (handles both single and multi-sequence structures)
        annotation_length = gene_info.get('annotated_gene_length', 'Unknown')
        copies = gene_info['copies']
        median_length = gene_info.get('ref_median_length')
        translation_details = gene_info.get('translation_details', '')
        
        # Use existing warning data instead of recalculating
        issue = gene_info.get('issue', 'unknown')
        details = gene_info.get('details', '')
        translation_status = gene_info.get('translation_status', 'OK')
        
        # Build warning message from existing data
        warning_parts = []
        if issue == 'too_short' or issue == 'too_long' or issue == 'not_in_median_lengths':
            warning_parts.append(details)
        elif issue == 'OK':
            pass
        else:
            warning_parts.append('Unknown issue')
        
        if translation_status == 'FAIL':
            warning_parts.append('Translation issue')
        
        warning_msg = ', '.join(warning_parts) if warning_parts else 'OK'
        
        return {
            'sequence_name': sequence_name,
            'gene_name': gene_name,
            'annotation_length': annotation_length,
            'median_length': median_length,
            'copies': copies,
            'warning_msg': warning_msg,
            'translation_details': translation_details
        }
    
    def add_missing_genes_to_reports(missing_genes, sequence_names, sample_tsv_lines, combined_tsv_lines):
        """Add missing genes to both individual and combined reports."""
        missing_count = 0
        for missing_gene in missing_genes:
            median_length = gene_median_lengths.get(missing_gene, 'NA')
            for seq_name in sequence_names:
                missing_row = format_tsv_row(
                    seq_name, missing_gene, 'NA', median_length, min_threshold, 
                    max_threshold, 'NA', 'Missing', 'NA'
                )
                combined_tsv_lines.append('\t'.join(missing_row))
                sample_tsv_lines.append('\t'.join(missing_row))
                missing_count += 1
        return missing_count

    logger.info(f"{"[INFO]:":10} Writing gene validation reports to: {report_directory}\n")
    
    # Prepare combined TSV data
    combined_tsv_lines = ['\t'.join(create_tsv_header())]
    total_warnings_in_reports = 0

    # Process each sample
    for sample_name, result in all_results.items():
        if 'error' in result:
            # Add error row to combined report
            error_row = format_tsv_row(
                sample_name, 'ERROR', 0, 0, min_threshold, max_threshold, 
                0, result['error'], 'ERROR'
            )
            combined_tsv_lines.append('\t'.join(error_row))
            continue
        
        # Get sequence data
        sequence_gene_data = result.get('sequence_gene_data', {})
        
        # Count total genes across all sequences
        total_genes_count = sum(len(seq_data['sequence_gene_data']) for seq_data in sequence_gene_data.values())
        logger.debug(f"{"[DEBUG]:":10} Processing {sample_name}: found {total_genes_count} genes across "
                     f"{len(sequence_gene_data)} sequences")
        
        # Prepare individual sample report
        sample_tsv_lines = ['\t'.join(create_tsv_header())]
        
        warnings_found = 0
        missing_genes_count = 0
        
        # Process all genes from sequence data
        sequence_names = list(sequence_gene_data.keys())
        for seq_name, seq_data in sequence_gene_data.items():
            gene_info_dict = seq_data['sequence_gene_data']
            
            for gene_name, gene_info in sorted(gene_info_dict.items()):
                gene_data = process_gene_data(gene_name, gene_info, seq_name)
                
                # Count warnings
                if gene_data['warning_msg'] != 'OK':
                    warnings_found += 1
                    total_warnings_in_reports += 1
                
                # Add to reports
                row = format_tsv_row(
                    gene_data['sequence_name'], gene_data['gene_name'], 
                    gene_data['annotation_length'], gene_data['median_length'],
                    min_threshold, max_threshold, gene_data['copies'],
                    gene_data['warning_msg'],
                    gene_data['translation_details']
                )
                combined_tsv_lines.append('\t'.join(row))
                sample_tsv_lines.append('\t'.join(row))
        
        # Add missing genes
        missing_genes = result.get('missing_genes', [])
        missing_genes_count = add_missing_genes_to_reports(
            missing_genes, sequence_names, sample_tsv_lines, combined_tsv_lines
        )
        total_warnings_in_reports += missing_genes_count
        
        # Write individual sample report
        sample_report_file = os.path.join(report_directory, f"{sample_name}_gene_validation_report.tsv")
        with open(sample_report_file, 'w') as f:
            f.write('\n'.join(sample_tsv_lines))
        
        logger.info(f"{"":10} {sample_name:20}: {total_genes_count} genes (including multi-copy), "
                    f"{warnings_found} warnings, {missing_genes_count} missing genes (from 113-gene reference set)")
    
    logger.info(f"")
    
    # Write combined TSV file
    combined_report_file = os.path.join(report_directory, "all_samples_gene_validation_report.tsv")
    with open(combined_report_file, 'w') as f:
        f.write('\n'.join(combined_tsv_lines))


def get_references(args, gene_median_lengths, gene_synonyms=None, data_dir_base=None):
    """
    Check alignment settings and get reference sequences from appropriate sources.
    
    Args:
        args: Command line arguments containing refs_order, refs_folder, and no_alignment
        gene_median_lengths (dict): Dictionary of median gene lengths
        gene_synonyms (dict): Dictionary mapping gene names to standardized synonyms
        
    Returns:
        dict: Dictionary with gene names as keys and lists of SeqRecord objects as values
    """
    
    # Initialize reference sources
    order_refs = None
    default_refs = None
    custom_refs = None
    
    # Require provided data_dir_base and do not fall back
    if data_dir_base is None:
        raise ValueError("data_dir_base must be provided to get_references()")

    # Check if refs_order contains entries
    if args.refs_order:
        # If refs_order is specified, no_alignment must be False
        if args.no_alignment:
            logger.error(f"{"[ERROR]:":10} Cannot specify --refs_order when --no_alignment is True. Please remove "
                         f"--no_alignment or remove --refs_order.")
            utils.exit_program()
        
        # Get available order directories under the resolved base data directory
        order_genomes_dir = os.path.join(data_dir_base, 'order_genomes')
        
        if not os.path.exists(order_genomes_dir):
            logger.error(f"Order genomes directory not found: {order_genomes_dir}")
            utils.exit_program()
        
        # Get list of available order directories
        available_orders = [d.lower() for d in os.listdir(order_genomes_dir) 
                          if os.path.isdir(os.path.join(order_genomes_dir, d))]
        
        # Check if all requested orders are available
        missing_orders = []
        for order in args.refs_order:
            if order.lower() not in available_orders:
                missing_orders.append(order)
        
        if missing_orders:
            logger.error(f"Reference plastid genomes not found for the following orders: {missing_orders}")
            logger.error(f"Available orders: {sorted(available_orders)}")
            utils.exit_program()
        
        # Get reference sequences from order-specific directories
        logger.info(f"{"[INFO]:":10} Using reference genomes from orders: {args.refs_order}")
        order_refs = get_ref_gene_seqrecords_from_orders(
            args.refs_order, gene_median_lengths, gene_synonyms, data_dir_base=data_dir_base
        )
    
    # Check if custom reference folder is specified
    if args.custom_refs_folder:
        # If custom folder is specified, no_alignment must be False
        if args.no_alignment:
            logger.error(f"{"[ERROR]:":10} Cannot specify --custom_refs_folder when --no_alignment is True. Please "
                         f"remove --no_alignment or remove --custom_refs_folder.")
            utils.exit_program()
        
        try:
            logger.info(f"{"[INFO]:":10} Using custom reference genomes from: {args.custom_refs_folder}")
            custom_refs = get_ref_gene_seqrecords_from_custom_folder(args.custom_refs_folder, gene_median_lengths,
                                                                     gene_synonyms)
        except Exception as e:
            logger.error(f"{"[ERROR]:":10} Error processing custom reference folder: {e}")
            utils.exit_program()
    
    # If no specific references are specified, use default reference genomes
    if not args.refs_order and not args.custom_refs_folder:
        logger.info(f"{"[INFO]:":10} Using default reference genomes from data/reference_genomes_default")
        default_refs = get_ref_gene_seqrecords_from_default(gene_median_lengths, gene_synonyms, 
                                                            data_dir_base=data_dir_base)
    
    # Merge all reference sources
    merged_refs = merge_reference_sequences(order_refs, default_refs, custom_refs)
    
    # Log summary of reference sources used
    sources_used = []
    if order_refs:
        sources_used.append("order-specific references")
    if default_refs:
        sources_used.append("default references")
    if custom_refs:
        sources_used.append("custom references")
    
    if sources_used:
        logger.info(f"{"[INFO]:":10} Combined reference sequences from: {', '.join(sources_used)}")

    utils.log_separator(logger)
    
    return merged_refs


def extract_taxonomic_info(record, order_name=None):
    """
    Extract taxonomic information from GenBank record for ID generation.
    
    Args:
        record: BioPython SeqRecord object
        order_name: The order name from the directory structure
        
    Returns:
        tuple: (order, family, genus, species) for ID generation
    """
    # Get organism name
    organism = record.annotations.get('organism', 'Unknown')
    
    # Extract genus and species from organism name
    organism_parts = organism.split()
    if len(organism_parts) >= 2:
        genus = organism_parts[0]
        species = organism_parts[1]
    else:
        genus = organism
        species = 'sp'
    
    # Get family from taxonomy
    taxonomy = record.annotations.get('taxonomy', [])
    family = 'Unknown'
    
    # Look for family in taxonomy (usually ends with 'aceae')
    for taxon in taxonomy:
        if taxon.endswith('aceae'):
            family = taxon
            break
    
    # Get order from taxonomy or use provided order name
    order = 'Unknown'
    if order_name:
        order = order_name
    else:
        # Look for order in taxonomy (usually ends with 'ales')
        for taxon in taxonomy:
            if taxon.endswith('ales'):
                order = taxon
                break
    
    return order, family, genus, species


def get_ref_gene_seqrecords_from_orders(orders, gene_median_lengths, gene_synonyms=None, data_dir_base=None):
    """
    Extract CDS, rRNA, and tRNA sequences from GenBank files in order-specific directories.
    
    Args:
        orders (list): List of order names to use for references
        gene_median_lengths (dict): Dictionary of median gene lengths
        gene_synonyms (dict): Dictionary mapping gene names to standardized synonyms
        
    Returns:
        dict: Dictionary with gene names as keys and lists of SeqRecord objects as values.
        Multiple copies of the same gene from a single reference file are labeled with
        _copy_1, _copy_2, etc. in their sequence IDs.
    """
    ref_cds_seqrecords = {}
    
    gene_names = set()
    mapped_genes = set()
    unmapped_genes = set()

    chloe_gene_names = gene_median_lengths.keys()
    
    # Initialize synonyms if not provided
    if gene_synonyms is None:
        gene_synonyms = {}

    # Get path to the order genomes data directory using base data dir 
    order_genomes_dir = os.path.join(data_dir_base, 'order_genomes')
    
    # Track gene copies per file for ID generation
    gene_copy_counts = {}  # (mapped_gene_name, gbk_file) -> copy_count
    
    for order in orders:
        order_dir = os.path.join(order_genomes_dir, order)
        gbk_files = glob.glob(os.path.join(order_dir, '*.gb*'))
        
        logger.debug(f"{"[DEBUG]:":10} Found {len(gbk_files)} GenBank files in {order} directory")
        
        for gbk_file in gbk_files:
            try:
                records = utils.parse_genbank_file(gbk_file, logger)
                
                for record in records:
                    for feature in record.features:
                        if feature.type in ["CDS", "rRNA", "tRNA"]:
                            gene_name = None
                            
                            # Extract gene name from qualifiers
                            if 'gene' in feature.qualifiers:
                                gene_name = feature.qualifiers['gene'][0]
                            elif 'product' in feature.qualifiers:
                                gene_name = '_'.join(feature.qualifiers['product'])
                            else:
                                raise ValueError(f"No gene name found for feature: {feature}")
                            
                            if gene_name:
                                gene_names.add(gene_name)
                                
                                # Try to map gene name using synonyms
                                mapped_gene_name = gene_name
                                if gene_name in gene_synonyms:
                                    mapped_gene_name = gene_synonyms[gene_name]
                                    mapped_genes.add(gene_name)
                                    logger.debug(f"{"[DEBUG]:":10} Mapped gene name: {gene_name} -> {mapped_gene_name}")
                                else:
                                    unmapped_genes.add(gene_name)
                                
                                # Extract the gene sequence
                                gene_seq = feature.location.extract(record.seq)
                                
                                # Extract taxonomic information for ID generation
                                order_tax, family_tax, genus_tax, species_tax = extract_taxonomic_info(record, order)
                                
                                # Track copy count for this gene in this file
                                copy_key = (mapped_gene_name, gbk_file)
                                if copy_key not in gene_copy_counts:
                                    gene_copy_counts[copy_key] = 0
                                gene_copy_counts[copy_key] += 1
                                copy_number = gene_copy_counts[copy_key]
                                
                                # Create a SeqRecord for this gene with copy number
                                base_id = f"{order_tax}_{family_tax}_{genus_tax}_{species_tax}_{mapped_gene_name}_{os.path.basename(gbk_file)}"
                                if copy_number > 1:
                                    gene_id = f"{base_id}_copy_{copy_number}"
                                else:
                                    gene_id = base_id
                                
                                gene_seqrecord = SeqRecord(
                                    seq=gene_seq,
                                    id=gene_id,
                                    description=f"{feature.type} from {order}/{os.path.basename(gbk_file)} (original: {gene_name})"
                                )
                                
                                # Add to the dictionary using mapped name
                                if mapped_gene_name not in ref_cds_seqrecords:
                                    ref_cds_seqrecords[mapped_gene_name] = []
                                ref_cds_seqrecords[mapped_gene_name].append(gene_seqrecord)
                
            except Exception as e:
                logger.warning(f"{"[WARNING]:":10} Error processing GenBank file {os.path.basename(gbk_file)}: {e}")
                continue
    
    # Check which Genbank gene names are still not found in the Chloe gene names after mapping
    unfound_gbk_gene_names = []
    for gk_gene_name in gene_names:
        mapped_name = gene_synonyms.get(gk_gene_name, gk_gene_name)
        if mapped_name not in chloe_gene_names:
            unfound_gbk_gene_names.append(gk_gene_name)

    if unfound_gbk_gene_names:
        logger.debug(f"{"[DEBUG]:":10} Genes still not found in median lengths after mapping: {unfound_gbk_gene_names}...")
    
    # Log summary
    total_genes = len(ref_cds_seqrecords)
    total_sequences = sum(len(seqrecords) for seqrecords in ref_cds_seqrecords.values())
    logger.debug(f"{"[DEBUG]:":10} Extracted {total_sequences} gene sequences for {total_genes} unique genes from "
                 f"order-specific GenBank files")
    
    return ref_cds_seqrecords


def merge_reference_sequences(*reference_dicts):
    """
    Merge multiple reference sequence dictionaries into a single dictionary.
    
    Args:
        *reference_dicts: Variable number of reference sequence dictionaries
        
    Returns:
        dict: Merged dictionary with gene names as keys and concatenated lists of SeqRecord objects as values
    """
    merged_refs = {}
    
    for ref_dict in reference_dicts:
        if ref_dict is None:
            continue
            
        for gene_name, seqrecords in ref_dict.items():
            if gene_name not in merged_refs:
                merged_refs[gene_name] = []
            merged_refs[gene_name].extend(seqrecords)
    
    return merged_refs


def get_ref_gene_seqrecords_from_custom_folder(refs_folder, gene_median_lengths, gene_synonyms=None):
    """
    Extract CDS, rRNA, and tRNA sequences from GenBank files in a custom reference directory.
    
    Args:
        refs_folder (str): Path to custom reference folder containing GenBank files
        gene_median_lengths (dict): Dictionary of median gene lengths
        gene_synonyms (dict): Dictionary mapping gene names to standardized synonyms
        
    Returns:
        dict: Dictionary with gene names as keys and lists of SeqRecord objects as values.
        Multiple copies of the same gene from a single reference file are labeled with
        _copy_1, _copy_2, etc. in their sequence IDs.
    """
    ref_cds_seqrecords = {}
    
    gene_names = set()
    mapped_genes = set()
    unmapped_genes = set()

    chloe_gene_names = gene_median_lengths.keys()
    
    # Initialize synonyms if not provided
    if gene_synonyms is None:
        gene_synonyms = {}

    # Check if the custom folder exists
    if not os.path.exists(refs_folder):
        raise ValueError(f"Custom reference folder not found: {refs_folder}")
    
    if not os.path.isdir(refs_folder):
        raise ValueError(f"Custom reference path is not a directory: {refs_folder}")
    
    # Get GenBank files from the custom folder
    gbk_files = glob.glob(os.path.join(refs_folder, '*.gb*'))
    
    logger.debug(f"{"[DEBUG]:":10} Found {len(gbk_files)} GenBank files in custom reference folder: {refs_folder}")
    
    # Track gene copies per file for ID generation
    gene_copy_counts = {}  # (mapped_gene_name, gbk_file) -> copy_count
    
    for gbk_file in gbk_files:
        try:
            records = utils.parse_genbank_file(gbk_file, logger)
            
            for record in records:
                for feature in record.features:
                    if feature.type in ["CDS", "rRNA", "tRNA"]:
                        gene_name = None
                        
                        # Extract gene name from qualifiers
                        if 'gene' in feature.qualifiers:
                            gene_name = feature.qualifiers['gene'][0]
                        elif 'product' in feature.qualifiers:
                            gene_name = '_'.join(feature.qualifiers['product'])
                        else:
                            raise ValueError(f"No gene name found for feature: {feature}")
                        
                        if gene_name:
                            gene_names.add(gene_name)
                            
                            # Try to map gene name using synonyms
                            mapped_gene_name = gene_name
                            if gene_name in gene_synonyms:
                                mapped_gene_name = gene_synonyms[gene_name]
                                mapped_genes.add(gene_name)
                                logger.debug(f"{"[DEBUG]:":10} Mapped gene name: {gene_name} -> {mapped_gene_name}")
                            else:
                                unmapped_genes.add(gene_name)
                            
                            # Extract the gene sequence
                            gene_seq = feature.location.extract(record.seq)
                            
                            # Extract taxonomic information for ID generation
                            order_tax, family_tax, genus_tax, species_tax = extract_taxonomic_info(record)
                            
                            # Track copy count for this gene in this file
                            copy_key = (mapped_gene_name, gbk_file)
                            if copy_key not in gene_copy_counts:
                                gene_copy_counts[copy_key] = 0
                            gene_copy_counts[copy_key] += 1
                            copy_number = gene_copy_counts[copy_key]
                            
                            # Create a SeqRecord for this gene with copy number
                            base_id = f"{order_tax}_{family_tax}_{genus_tax}_{species_tax}_{mapped_gene_name}_{os.path.basename(gbk_file)}"
                            if copy_number > 1:
                                gene_id = f"{base_id}_copy_{copy_number}"
                            else:
                                gene_id = base_id
                            
                            gene_seqrecord = SeqRecord(
                                seq=gene_seq,
                                id=gene_id,
                                description=f"{feature.type} from custom/{os.path.basename(gbk_file)} (original: {gene_name})"
                            )
                            
                            # Add to the dictionary using mapped name
                            if mapped_gene_name not in ref_cds_seqrecords:
                                ref_cds_seqrecords[mapped_gene_name] = []
                            ref_cds_seqrecords[mapped_gene_name].append(gene_seqrecord)
            
        except Exception as e:
            logger.warning(f"{"[WARNING]:":10} Error processing GenBank file {os.path.basename(gbk_file)}: {e}")
            continue
    
    # Check which Genbank gene names are still not found in the Chloe gene names after mapping
    unfound_gbk_gene_names = []
    for gk_gene_name in gene_names:
        mapped_name = gene_synonyms.get(gk_gene_name, gk_gene_name)
        if mapped_name not in chloe_gene_names:
            unfound_gbk_gene_names.append(gk_gene_name)

    if unfound_gbk_gene_names:
        logger.debug(f"{"[DEBUG]:":10} Genes still not found in median lengths after mapping: {unfound_gbk_gene_names}...")
    
    # Log summary
    total_genes = len(ref_cds_seqrecords)
    total_sequences = sum(len(seqrecords) for seqrecords in ref_cds_seqrecords.values())
    logger.debug(f"{"[DEBUG]:":10} Extracted {total_sequences} gene sequences for {total_genes} unique genes from custom reference folder")
    
    return ref_cds_seqrecords


def get_ref_gene_seqrecords_from_default(gene_median_lengths, gene_synonyms=None, data_dir_base=None):
    """
    Extract CDS, rRNA, and tRNA sequences from GenBank files in the default reference directory.
    
    Args:
        gene_median_lengths (dict): Dictionary of median gene lengths
        gene_synonyms (dict): Dictionary mapping gene names to standardized synonyms
        
    Returns:
        dict: Dictionary with gene names as keys and lists of SeqRecord objects as values.
        Multiple copies of the same gene from a single reference file are labeled with
        _copy_1, _copy_2, etc. in their sequence IDs.
    """
    ref_cds_seqrecords = {}
    
    gene_names = set()
    mapped_genes = set()
    unmapped_genes = set()

    chloe_gene_names = gene_median_lengths.keys()
    
    # Initialize synonyms if not provided
    if gene_synonyms is None:
        gene_synonyms = {}

    # Get path to the default reference genomes directory using base data dir
    ref_geneome_dir = os.path.join(data_dir_base, 'reference_genomes_default')
    
    # Track gene copies per file for ID generation
    gene_copy_counts = {}  # (mapped_gene_name, gbk_file) -> copy_count
    gbk_files = glob.glob(os.path.join(ref_geneome_dir, '*.gb*'))
    
    logger.debug(f"{"[DEBUG]:":10} Found {len(gbk_files)} GenBank files in default reference data directory")
    
    for gbk_file in gbk_files:
        try:
            records = utils.parse_genbank_file(gbk_file, logger)
            
            for record in records:
                for feature in record.features:
                    if feature.type in ["CDS", "rRNA", "tRNA"]:
                        gene_name = None
                        
                        # Extract gene name from qualifiers
                        if 'gene' in feature.qualifiers:
                            gene_name = feature.qualifiers['gene'][0]
                        elif 'product' in feature.qualifiers:
                            gene_name = '_'.join(feature.qualifiers['product'])
                        else:
                            raise ValueError(f"No gene name found for feature: {feature}")
                        
                        if gene_name:
                            gene_names.add(gene_name)
                            
                            # Try to map gene name using synonyms
                            mapped_gene_name = gene_name
                            if gene_name in gene_synonyms:
                                mapped_gene_name = gene_synonyms[gene_name]
                                mapped_genes.add(gene_name)
                                logger.debug(f"{"[DEBUG]:":10} Mapped gene name: {gene_name} -> {mapped_gene_name}")
                            else:
                                unmapped_genes.add(gene_name)
                            
                            # Extract the gene sequence
                            gene_seq = feature.location.extract(record.seq)
                            
                            # Extract taxonomic information for ID generation
                            order_tax, family_tax, genus_tax, species_tax = extract_taxonomic_info(record)
                            
                            # Track copy count for this gene in this file
                            copy_key = (mapped_gene_name, gbk_file)
                            if copy_key not in gene_copy_counts:
                                gene_copy_counts[copy_key] = 0
                            gene_copy_counts[copy_key] += 1
                            copy_number = gene_copy_counts[copy_key]
                            
                            # Create a SeqRecord for this gene with copy number
                            base_id = f"{order_tax}_{family_tax}_{genus_tax}_{species_tax}_{mapped_gene_name}_{os.path.basename(gbk_file)}"
                            if copy_number > 1:
                                gene_id = f"{base_id}_copy_{copy_number}"
                            else:
                                gene_id = base_id
                            
                            gene_seqrecord = SeqRecord(
                                seq=gene_seq,
                                id=gene_id,
                                description=f"{feature.type} from default/{os.path.basename(gbk_file)} (original: {gene_name})"
                            )
                            
                            # Add to the dictionary using mapped name
                            if mapped_gene_name not in ref_cds_seqrecords:
                                ref_cds_seqrecords[mapped_gene_name] = []
                            ref_cds_seqrecords[mapped_gene_name].append(gene_seqrecord)
            
        except Exception as e:
            logger.warning(f"{"[WARNING]:":10} Error processing GenBank file {os.path.basename(gbk_file)}: {e}")
            continue
    
    # Check which Genbank gene names are still not found in the Chloe gene names after mapping
    unfound_gbk_gene_names = []
    for gk_gene_name in gene_names:
        mapped_name = gene_synonyms.get(gk_gene_name, gk_gene_name)
        if mapped_name not in chloe_gene_names:
            unfound_gbk_gene_names.append(gk_gene_name)

    if unfound_gbk_gene_names:
        logger.debug(f"{"[DEBUG]:":10} Genes still not found in median lengths after mapping: {unfound_gbk_gene_names}...")
    
    # Log summary
    total_genes = len(ref_cds_seqrecords)
    total_sequences = sum(len(seqrecords) for seqrecords in ref_cds_seqrecords.values())
    logger.debug(f"{"[DEBUG]:":10} Extracted {total_sequences} gene sequences for {total_genes} unique genes from default reference GenBank files")
    
    return ref_cds_seqrecords


def create_single_cds_alignment(task_data):
    """
    Create alignment for a single gene from a single sample.
    
    Args:
        task_data (dict): Dictionary containing all necessary data for alignment
        
    Returns:
        tuple: (success, result) where success is bool and result is either success message or error info
    """
    try:
        sample_name = task_data['sample_name']
        gene_name = task_data['gene_name']
        cds_info_list = task_data['cds_info_list']
        ref_sequences = task_data['ref_sequences']
        outdir_alignments = task_data['outdir_alignments']
        threads = task_data['threads']

        # Check if alignment file already exists and is not empty
        alignment_file = os.path.join(outdir_alignments, f"{sample_name}_{gene_name}_alignment.fasta")
        if utils.file_exists_and_not_empty(alignment_file):
            return True, f"Skipped {gene_name} - alignment already exists"

        seqrecords_to_write = []
    
        # Check for internal stop codons in sample sequences
        has_internal_stops = False

        for i, cds_info in enumerate(cds_info_list):
            seq = cds_info['cds_seq']
            # Convert sequence to Seqrecord for translation
            seqrecord = SeqRecord(seq=seq, id=f"{sample_name}_{gene_name}_copy_{i+1}", description=f"CDS from {sample_name}")
            seqrecords_to_write.append(seqrecord)
            translated = utils.pad_seq(seqrecord)[0].translate()
            if "*" in str(translated.seq)[:-1]:  # Exclude the last position
                has_internal_stops = True

        # Add reference sequences to sample seqrecords
        seqrecords_to_write.extend(ref_sequences)

        if has_internal_stops:
            # Align nucleotide sequences directly using mafft
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_nucleotide:

                SeqIO.write(seqrecords_to_write, temp_nucleotide, "fasta")
                temp_nucleotide_path = temp_nucleotide.name
            
            try:
                # Run MAFFT alignment
                subprocess.run(['mafft', '--auto', '--thread', str(threads), '--anysymbol', temp_nucleotide_path], 
                             stdout=open(alignment_file, 'w'), 
                             stderr=subprocess.PIPE, check=True)
                
                return True, f"Created nucleotide alignment for {gene_name}"
                
            except Exception as e:
                return False, (e, traceback.format_exc())
            
            finally:
                if os.path.exists(temp_nucleotide_path):
                    os.unlink(temp_nucleotide_path)
                
        else:
            # Write ungapped nucleotide seqs for backtranslation
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_nucleotide:
                for seqrecord in seqrecords_to_write:
                    seqrecord.seq = seqrecord.seq.replace("-", "")

                SeqIO.write(seqrecords_to_write, temp_nucleotide, "fasta")
                temp_nucleotide_path = temp_nucleotide.name

            # Create protein alignment using MAFFT, then backtranslate
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_protein:
                # Translate all seqrecords to protein
                for seqrecord in seqrecords_to_write:
                    seqrecord = utils.pad_seq(seqrecord)[0]
                    seqrecord.seq = seqrecord.seq.translate()

                SeqIO.write(seqrecords_to_write, temp_protein, "fasta")
                temp_protein_path = temp_protein.name

            try:
                # Align proteins using MAFFT
                aligned_protein_file = os.path.join(outdir_alignments, f"{sample_name}_{gene_name}_protein_aligned.fasta")
                subprocess.run(['mafft', '--auto', '--thread', str(threads), '--anysymbol', temp_protein_path], 
                             stdout=open(aligned_protein_file, 'w'), 
                             stderr=subprocess.PIPE, check=True)
                
                # Backtranslate using trimal
                subprocess.run(['trimal', '-backtrans', temp_nucleotide_path, '-in', aligned_protein_file, '-out', alignment_file, '-ignorestopcodon'], 
                             stderr=subprocess.PIPE, check=True)
                
                return True, f"Created backtranslated alignment for {gene_name}"
                
            except Exception as e:
                return False, (e, traceback.format_exc())
            
            finally:
                if os.path.exists(temp_nucleotide_path):
                    os.unlink(temp_nucleotide_path)
                if os.path.exists(temp_protein_path):
                    os.unlink(temp_protein_path)
                if os.path.exists(aligned_protein_file):
                    os.unlink(aligned_protein_file)
        
    except Exception as e:
        return False, (e, traceback.format_exc())


def generate_cds_alignments(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads):
    """
    Create alignments between sample CDS sequences and reference CDS sequences using multiprocessing.
    
    Args:
        all_sample_results (dict): Dictionary of all sample results
        ref_gene_seqrecords (dict): Dictionary of reference gene sequences
        logger: Logger instance for logging messages
        output_directory (str): Directory to write alignment files
        pool_size (int): Number of processes to use for multiprocessing
        threads (int): Number of threads to use for each process

    Returns:
        None
    """
    logger.debug(f"{"[DEBUG]:":10} Starting CDS alignment generation")

    
    # Create alignments directory with subfolder structure
    outdir_alignments = os.path.join(output_directory, '03_alignments_with_refs', '01_per_sample_alignments')
    os.makedirs(outdir_alignments, exist_ok=True)
    
    # Prepare alignment tasks
    alignment_tasks = []
    for sample_name, result in all_sample_results.items():
        if 'error' in result:
            continue
            
        # Create sample-specific subfolder
        sample_alignments_dir = os.path.join(outdir_alignments, sample_name)
        os.makedirs(sample_alignments_dir, exist_ok=True)
            
        gene_cds_info = result.get('gene_cds_info', {})
        
        # Group genes by base name (without copy suffix) to handle multiple copies
        gene_groups = {}
        for gene_name, cds_info_list in gene_cds_info.items():
            # Extract base gene name (remove _copy_N suffix if present)
            base_gene_name = gene_name.split('_copy_')[0]
            
            if base_gene_name not in gene_groups:
                gene_groups[base_gene_name] = []
            gene_groups[base_gene_name].extend(cds_info_list)
        
        for base_gene_name, all_cds_info in gene_groups.items():
            # Only create alignments for genes that are in reference
            if base_gene_name in ref_gene_seqrecords:
                ref_sequences = ref_gene_seqrecords[base_gene_name]

                alignment_tasks.append({
                    'sample_name': sample_name,
                    'gene_name': base_gene_name,
                    'cds_info_list': all_cds_info,
                    'ref_sequences': ref_sequences,
                    'outdir_alignments': sample_alignments_dir,  # Use sample-specific directory
                    'threads': threads
                })
            else:
                logger.warning(f"{"[WARNING]:":10} No reference sequences found for gene {base_gene_name}")
    
    if not alignment_tasks:
        logger.info(f"{"[INFO]:":10} No CDS alignment tasks found")
        return
    
    logger.debug(f"{"[DEBUG]:":10} Found {len(alignment_tasks)} CDS alignment tasks")
    
    # Process alignments with multiprocessing and progress bar
    logger.info(f"{"[INFO]:":10} Processing {len(alignment_tasks)} CDS alignment tasks with {pool_size} processes and {threads} threads per process")
    with ProcessPoolExecutor(max_workers=pool_size) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(create_single_cds_alignment, task): f"{task['sample_name']}_{task['gene_name']}"
            for task in alignment_tasks
        }
        
        # Process results with progress bar
        for future in tqdm(as_completed(future_to_task), total=len(future_to_task), desc=f"{"[INFO]:":10} {"Generating per-sample CDS alignments":<40}", file=sys.stdout):
            
            task_id = future_to_task[future]

            success, result = future.result()
            
            if success:
                logger.debug(f"{"[DEBUG]:":10} {task_id}: {result}")
            else:
                utils.log_manager.handle_error(result[0], result[1], "create_single_cds_alignment()", task_id)

    logger.debug(f"{"[DEBUG]:":10} CDS alignment generation complete")


def create_single_rRNA_alignment(task_data):
    """
    Create alignment for a single rRNA gene from a single sample.
    
    Args:
        task_data (dict): Dictionary containing all necessary data for alignment
        
    Returns:
        tuple: (success, result) where success is bool and result is either success message or error info
    """
    try:
        sample_name = task_data['sample_name']
        gene_name = task_data['gene_name']
        rRNA_info_list = task_data['rRNA_info_list']
        ref_sequences = task_data['ref_sequences']
        outdir_alignments = task_data['outdir_alignments']
        threads = task_data['threads']

        # Check if alignment file already exists and is not empty
        alignment_file = os.path.join(outdir_alignments, f"{sample_name}_{gene_name}_rRNA_alignment.fasta")
        if utils.file_exists_and_not_empty(alignment_file):
            return True, f"Skipped {gene_name} - rRNA alignment already exists"

        seqrecords_to_write = []
    
        # Add sample rRNA sequences
        for i, rRNA_info in enumerate(rRNA_info_list):
            seq = rRNA_info['rRNA_seq']
            # Convert sequence to Seqrecord
            seqrecord = SeqRecord(seq=seq, id=f"{sample_name}_{gene_name}_copy_{i+1}", description=f"rRNA from {sample_name}")
            seqrecords_to_write.append(seqrecord)

        # Add reference sequences
        seqrecords_to_write.extend(ref_sequences)

        # Align nucleotide sequences directly using mafft
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_nucleotide:
            SeqIO.write(seqrecords_to_write, temp_nucleotide, "fasta")
            temp_nucleotide_path = temp_nucleotide.name
        
        try:
            # Run MAFFT alignment
            subprocess.run(['mafft', '--auto', '--thread', str(threads), '--anysymbol', temp_nucleotide_path], 
                         stdout=open(alignment_file, 'w'), 
                         stderr=subprocess.PIPE, check=True)
            
            return True, f"Created rRNA alignment for {gene_name}"
            
        except Exception as e:
            return False, (e, traceback.format_exc())
        
        finally:
            if os.path.exists(temp_nucleotide_path):
                os.unlink(temp_nucleotide_path)
        
    except Exception as e:
        return False, (e, traceback.format_exc())


def create_single_tRNA_alignment(task_data):
    """
    Create alignment for a single tRNA gene from a single sample.
    
    Args:
        task_data (dict): Dictionary containing all necessary data for alignment
        
    Returns:
        tuple: (success, result) where success is bool and result is either success message or error info
    """
    try:
        sample_name = task_data['sample_name']
        gene_name = task_data['gene_name']
        tRNA_info_list = task_data['tRNA_info_list']
        ref_sequences = task_data['ref_sequences']
        outdir_alignments = task_data['outdir_alignments']
        threads = task_data['threads']

        # Check if alignment file already exists and is not empty
        alignment_file = os.path.join(outdir_alignments, f"{sample_name}_{gene_name}_tRNA_alignment.fasta")
        if utils.file_exists_and_not_empty(alignment_file):
            return True, f"Skipped {gene_name} - tRNA alignment already exists"

        seqrecords_to_write = []
    
        # Add sample tRNA sequences
        for i, tRNA_info in enumerate(tRNA_info_list):
            seq = tRNA_info['tRNA_seq']
            # Convert sequence to Seqrecord
            seqrecord = SeqRecord(seq=seq, id=f"{sample_name}_{gene_name}_copy_{i+1}", description=f"tRNA from {sample_name}")
            seqrecords_to_write.append(seqrecord)

        # Add reference sequences
        seqrecords_to_write.extend(ref_sequences)

        # Align nucleotide sequences directly using mafft
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_nucleotide:
            SeqIO.write(seqrecords_to_write, temp_nucleotide, "fasta")
            temp_nucleotide_path = temp_nucleotide.name
        
        try:
            # Run MAFFT alignment
            subprocess.run(['mafft', '--auto', '--thread', str(threads), '--anysymbol', temp_nucleotide_path], 
                         stdout=open(alignment_file, 'w'), 
                         stderr=subprocess.PIPE, check=True)
            
            return True, f"Created tRNA alignment for {gene_name}"
            
        except Exception as e:
            return False, (e, traceback.format_exc())
        
        finally:
            if os.path.exists(temp_nucleotide_path):
                os.unlink(temp_nucleotide_path)
        
    except Exception as e:
        return False, (e, traceback.format_exc())


def generate_rRNA_alignments(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads):
    """
    Create alignments for rRNA genes using MAFFT for nucleotide alignments.
    
    Args:   
        all_sample_results (dict): Dictionary of all sample results
        ref_gene_seqrecords (dict): Dictionary of reference sequences (will be filtered for rRNA)
        output_directory (str): Directory to write alignment files
        pool_size (int): Number of processes to use for multiprocessing
        threads (int): Number of threads to use for each process

    Returns:
        None
    """
    logger.debug(f"{"[DEBUG]:":10} Starting rRNA alignment generation")
    
    # Create alignments directory with subfolder structure
    outdir_alignments = os.path.join(output_directory, '03_alignments_with_refs', '01_per_sample_alignments')
    os.makedirs(outdir_alignments, exist_ok=True)
    
    # Prepare alignment tasks
    alignment_tasks = []
    for sample_name, result in all_sample_results.items():
        if 'error' in result:
            continue
            
        # Create sample-specific subfolder
        sample_alignments_dir = os.path.join(outdir_alignments, sample_name)
        os.makedirs(sample_alignments_dir, exist_ok=True)
            
        gene_rRNA_info = result.get('gene_rRNA_info', {})
        
        # Group genes by base name (without copy suffix) to handle multiple copies
        gene_groups = {}
        for gene_name, rRNA_info_list in gene_rRNA_info.items():
            # Extract base gene name (remove _copy_N suffix if present)
            base_gene_name = gene_name.split('_copy_')[0]
            
            if base_gene_name not in gene_groups:
                gene_groups[base_gene_name] = []
            gene_groups[base_gene_name].extend(rRNA_info_list)
        
        for base_gene_name, all_rRNA_info in gene_groups.items():
            # Only create alignments for genes that are in reference
            if base_gene_name in ref_gene_seqrecords:
                ref_sequences = ref_gene_seqrecords[base_gene_name]

                alignment_tasks.append({
                    'sample_name': sample_name,
                    'gene_name': base_gene_name,
                    'rRNA_info_list': all_rRNA_info,
                    'ref_sequences': ref_sequences,
                    'outdir_alignments': sample_alignments_dir,  # Use sample-specific directory
                    'threads': threads
                })
            else:
                logger.debug(f"{"[DEBUG]:":10} No reference sequences found for rRNA gene {base_gene_name}")
    
    if not alignment_tasks:
        logger.warning(f"{"[WARNING]:":10} No rRNA alignment tasks found")
        return
    
    logger.debug(f"{"[DEBUG]:":10} Found {len(alignment_tasks)} rRNA alignment tasks")
    
    # Process alignments with multiprocessing and progress bar
    logger.info(f"{"[INFO]:":10} Processing {len(alignment_tasks)} rRNA alignment tasks with {pool_size} processes and {threads} threads per process")
    with ProcessPoolExecutor(max_workers=pool_size) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(create_single_rRNA_alignment, task): f"{task['sample_name']}_{task['gene_name']}"
            for task in alignment_tasks
        }
        
        # Process results with progress bar
        for future in tqdm(as_completed(future_to_task), total=len(future_to_task), desc=f"{"[INFO]:":10} {"Generating per-sample rRNA alignments":<40}", file=sys.stdout):
            
            task_id = future_to_task[future]

            success, result = future.result()
            
            if success:
                logger.debug(f"{"[DEBUG]:":10} {task_id}: {result}")
            else:
                utils.log_manager.handle_error(result[0], result[1], "create_single_rRNA_alignment()", task_id)

    logger.debug(f"{"[DEBUG]:":10} rRNA alignment generation complete")


def generate_tRNA_alignments(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads):
    """
    Create alignments for tRNA genes using MAFFT for nucleotide alignments.
    
    Args:   
        all_sample_results (dict): Dictionary of all sample results
        ref_gene_seqrecords (dict): Dictionary of reference sequences (will be filtered for tRNA)
        output_directory (str): Directory to write alignment files
        pool_size (int): Number of processes to use for multiprocessing
        threads (int): Number of threads to use for each process

    Returns:
        None
    """
    logger.debug(f"{"[DEBUG]:":10} Starting tRNA alignment generation")
    
    # Create alignments directory with subfolder structure
    outdir_alignments = os.path.join(output_directory, '03_alignments_with_refs', '01_per_sample_alignments')
    os.makedirs(outdir_alignments, exist_ok=True)
    
    # Prepare alignment tasks
    alignment_tasks = []
    for sample_name, result in all_sample_results.items():
        if 'error' in result:
            continue
            
        # Create sample-specific subfolder
        sample_alignments_dir = os.path.join(outdir_alignments, sample_name)
        os.makedirs(sample_alignments_dir, exist_ok=True)
            
        gene_tRNA_info = result.get('gene_tRNA_info', {})
        
        # Group genes by base name (without copy suffix) to handle multiple copies
        gene_groups = {}
        for gene_name, tRNA_info_list in gene_tRNA_info.items():
            # Extract base gene name (remove _copy_N suffix if present)
            base_gene_name = gene_name.split('_copy_')[0]
            
            if base_gene_name not in gene_groups:
                gene_groups[base_gene_name] = []
            gene_groups[base_gene_name].extend(tRNA_info_list)
        
        for base_gene_name, all_tRNA_info in gene_groups.items():
            # Only create alignments for genes that are in reference
            if base_gene_name in ref_gene_seqrecords:
                ref_sequences = ref_gene_seqrecords[base_gene_name]

                alignment_tasks.append({
                    'sample_name': sample_name,
                    'gene_name': base_gene_name,
                    'tRNA_info_list': all_tRNA_info,
                    'ref_sequences': ref_sequences,
                    'outdir_alignments': sample_alignments_dir,  # Use sample-specific directory
                    'threads': threads
                })
            else:
                logger.debug(f"{"[DEBUG]:":10} No reference sequences found for tRNA gene {base_gene_name}")
    
    if not alignment_tasks:
        logger.warning(f"{"[WARNING]:":10} No tRNA alignment tasks found")
        return
    
    logger.debug(f"{"[DEBUG]:":10} Found {len(alignment_tasks)} tRNA alignment tasks")
    
    # Process alignments with multiprocessing and progress bar
    logger.info(f"{"[INFO]:":10} Processing {len(alignment_tasks)} tRNA alignment tasks with {pool_size} processes and {threads} threads per process")
    with ProcessPoolExecutor(max_workers=pool_size) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(create_single_tRNA_alignment, task): f"{task['sample_name']}_{task['gene_name']}"
            for task in alignment_tasks
        }
        
        # Process results with progress bar
        for future in tqdm(as_completed(future_to_task), total=len(future_to_task), desc=f"{"[INFO]:":10} {"Generating per-sample tRNA alignments":<40}", file=sys.stdout):
            
            task_id = future_to_task[future]

            success, result = future.result()
            
            if success:
                logger.debug(f"{"[DEBUG]:":10} {task_id}: {result}")
            else:
                utils.log_manager.handle_error(result[0], result[1], "create_single_tRNA_alignment()", task_id)

    logger.debug(f"{"[DEBUG]:":10} tRNA alignment generation complete")


def create_single_per_gene_alignment(task_data):
    """
    Create alignment for a single gene across all samples.
    
    Args:
        task_data (dict): Dictionary containing all necessary data for alignment
        
    Returns:
        tuple: (success, result) where success is bool and result is either success message or error info
    """
    try:
        gene_name = task_data['gene_name']
        all_sample_sequences = task_data['all_sample_sequences']
        ref_sequences = task_data['ref_sequences']
        outdir_alignments = task_data['outdir_alignments']
        threads = task_data['threads']
        gene_type = task_data['gene_type']  # 'CDS', 'rRNA', or 'tRNA'

        # Check if alignment file already exists and is not empty
        alignment_file = os.path.join(outdir_alignments, f"{gene_name}_all_samples_alignment.fasta")
        if utils.file_exists_and_not_empty(alignment_file):
            return True, f"Skipped {gene_name} - per-gene alignment already exists"

        seqrecords_to_write = []
        
        # Add all sample sequences
        seqrecords_to_write.extend(all_sample_sequences)
        
        # Add reference sequences
        seqrecords_to_write.extend(ref_sequences)

        if gene_type == 'CDS':
            # Check for internal stop codons in sample sequences
            has_internal_stops = False
            for seqrecord in all_sample_sequences:
                translated = utils.pad_seq(seqrecord)[0].translate()
                if "*" in str(translated.seq)[:-1]:  # Exclude the last position
                    has_internal_stops = True
                    break
                    
            if has_internal_stops:
                # Align nucleotide sequences directly using mafft
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_nucleotide:
                    SeqIO.write(seqrecords_to_write, temp_nucleotide, "fasta")
                    temp_nucleotide_path = temp_nucleotide.name
                
                try:
                    # Run MAFFT alignment
                    subprocess.run(['mafft', '--auto', '--thread', str(threads), '--anysymbol', temp_nucleotide_path], 
                                 stdout=open(alignment_file, 'w'), 
                                 stderr=subprocess.PIPE, check=True)
                    
                    return True, f"Created nucleotide alignment for {gene_name} (all samples)"
                    
                except Exception as e:
                    return False, (e, traceback.format_exc())
                
                finally:
                    if os.path.exists(temp_nucleotide_path):
                        os.unlink(temp_nucleotide_path)
                    
            else:
                # Write ungapped nucleotide seqs for backtranslation
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_nucleotide:
                    for seqrecord in seqrecords_to_write:
                        seqrecord.seq = seqrecord.seq.replace("-", "")

                    SeqIO.write(seqrecords_to_write, temp_nucleotide, "fasta")
                    temp_nucleotide_path = temp_nucleotide.name

                # Create protein alignment using MAFFT, then backtranslate
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_protein:
                    # Translate all seqrecords to protein
                    for seqrecord in seqrecords_to_write:
                        seqrecord = utils.pad_seq(seqrecord)[0]
                        seqrecord.seq = seqrecord.seq.translate()

                    SeqIO.write(seqrecords_to_write, temp_protein, "fasta")
                    temp_protein_path = temp_protein.name

                try:
                    # Align proteins using MAFFT
                    aligned_protein_file = os.path.join(outdir_alignments, f"{gene_name}_protein_aligned.fasta")
                    subprocess.run(['mafft', '--auto', '--thread', str(threads), '--anysymbol', temp_protein_path], 
                                 stdout=open(aligned_protein_file, 'w'), 
                                 stderr=subprocess.PIPE, check=True)
                    
                    # Backtranslate using trimal
                    subprocess.run(['trimal', '-backtrans', temp_nucleotide_path, '-in', aligned_protein_file, '-out', alignment_file, '-ignorestopcodon'], 
                                 stderr=subprocess.PIPE, check=True)

                    return True, f"Created backtranslated alignment for {gene_name} (all samples)"
                    
                except Exception as e:
                    return False, (e, traceback.format_exc())
                
                finally:
                    if os.path.exists(temp_nucleotide_path):
                        os.unlink(temp_nucleotide_path)
                    if os.path.exists(temp_protein_path):
                        os.unlink(temp_protein_path)
                    if os.path.exists(aligned_protein_file):
                        os.unlink(aligned_protein_file)
        
        elif gene_type == 'rRNA':
            # Align nucleotide sequences directly using mafft
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_nucleotide:
                SeqIO.write(seqrecords_to_write, temp_nucleotide, "fasta")
                temp_nucleotide_path = temp_nucleotide.name
            
            try:
                # Run MAFFT alignment
                subprocess.run(['mafft', '--auto', '--thread', str(threads), '--anysymbol',temp_nucleotide_path], 
                             stdout=open(alignment_file, 'w'), 
                             stderr=subprocess.PIPE, check=True)
                
                return True, f"Created rRNA alignment for {gene_name} (all samples)"
                
            except Exception as e:
                return False, (e, traceback.format_exc())
            
            finally:
                if os.path.exists(temp_nucleotide_path):
                    os.unlink(temp_nucleotide_path)
        
        else:  # tRNA
            # Align nucleotide sequences directly using mafft
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_nucleotide:
                SeqIO.write(seqrecords_to_write, temp_nucleotide, "fasta")
                temp_nucleotide_path = temp_nucleotide.name
            
            try:
                # Run MAFFT alignment
                subprocess.run(['mafft', '--auto', '--thread', str(threads), '--anysymbol', temp_nucleotide_path], 
                             stdout=open(alignment_file, 'w'), 
                             stderr=subprocess.PIPE, check=True)
                
                return True, f"Created tRNA alignment for {gene_name} (all samples)"
                
            except Exception as e:
                return False, (e, traceback.format_exc())
            
            finally:
                if os.path.exists(temp_nucleotide_path):
                    os.unlink(temp_nucleotide_path)
        
    except Exception as e:
        return False, (e, traceback.format_exc())


def generate_per_gene_alignments(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads):
    """
    Create per-gene alignments that include all samples for each gene.
    
    Args:
        all_sample_results (dict): Dictionary of all sample results
        ref_gene_seqrecords (dict): Dictionary of reference gene sequences
        output_directory (str): Directory to write alignment files
        pool_size (int): Number of processes to use for multiprocessing
        threads (int): Number of threads to use for each process

    Returns:
        None
    """
    logger.debug(f"{"[DEBUG]:":10} Starting per-gene alignment generation")
    
    # Create alignments directory with subfolder structure
    outdir_alignments = os.path.join(output_directory, '03_alignments_with_refs', '02_per_gene_alignments')
    os.makedirs(outdir_alignments, exist_ok=True)
    
    # Collect all gene sequences from all samples
    gene_sequences = {}  # gene_name -> list of SeqRecord objects
    gene_types = {}      # gene_name -> 'CDS', 'rRNA', or 'tRNA'
    
    for sample_name, result in all_sample_results.items():
        if 'error' in result:
            continue
            
        gene_cds_info = result.get('gene_cds_info', {})
        gene_rRNA_info = result.get('gene_rRNA_info', {})
        gene_tRNA_info = result.get('gene_tRNA_info', {})
        
        # Process CDS genes
        for gene_name, cds_info_list in gene_cds_info.items():
            # Extract base gene name (remove _copy_N suffix if present)
            base_gene_name = gene_name.split('_copy_')[0]
            
            if base_gene_name not in gene_sequences:
                gene_sequences[base_gene_name] = []
                gene_types[base_gene_name] = 'CDS'
            
            for i, cds_info in enumerate(cds_info_list):
                seq = cds_info['cds_seq']
                # Extract original copy number from gene name if present
                if '_copy_' in gene_name:
                    copy_num = gene_name.split('_copy_')[1]
                else:
                    copy_num = str(i + 1)
                seqrecord = SeqRecord(seq=seq, id=f"{sample_name}_{base_gene_name}_copy_{copy_num}", description=f"CDS from {sample_name}")
                gene_sequences[base_gene_name].append(seqrecord)
        
        # Process rRNA genes
        for gene_name, rRNA_info_list in gene_rRNA_info.items():
            # Extract base gene name (remove _copy_N suffix if present)
            base_gene_name = gene_name.split('_copy_')[0]
            
            if base_gene_name not in gene_sequences:
                gene_sequences[base_gene_name] = []
                gene_types[base_gene_name] = 'rRNA'
            
            for i, rRNA_info in enumerate(rRNA_info_list):
                seq = rRNA_info['rRNA_seq']
                # Extract original copy number from gene name if present
                if '_copy_' in gene_name:
                    copy_num = gene_name.split('_copy_')[1]
                else:
                    copy_num = str(i + 1)
                seqrecord = SeqRecord(seq=seq, id=f"{sample_name}_{base_gene_name}_copy_{copy_num}", description=f"rRNA from {sample_name}")
                gene_sequences[base_gene_name].append(seqrecord)
        
        # Process tRNA genes
        for gene_name, tRNA_info_list in gene_tRNA_info.items():
            # Extract base gene name (remove _copy_N suffix if present)
            base_gene_name = gene_name.split('_copy_')[0]
            
            if base_gene_name not in gene_sequences:
                gene_sequences[base_gene_name] = []
                gene_types[base_gene_name] = 'tRNA'
            
            for i, tRNA_info in enumerate(tRNA_info_list):
                seq = tRNA_info['tRNA_seq']
                # Extract original copy number from gene name if present
                if '_copy_' in gene_name:
                    copy_num = gene_name.split('_copy_')[1]
                else:
                    copy_num = str(i + 1)
                seqrecord = SeqRecord(seq=seq, id=f"{sample_name}_{base_gene_name}_copy_{copy_num}", description=f"tRNA from {sample_name}")
                gene_sequences[base_gene_name].append(seqrecord)
    
    # Prepare alignment tasks
    alignment_tasks = []
    for gene_name, sample_sequences in gene_sequences.items():
        if gene_name in ref_gene_seqrecords:
            ref_sequences = ref_gene_seqrecords[gene_name]
            gene_type = gene_types[gene_name]

            alignment_tasks.append({
                'gene_name': gene_name,
                'all_sample_sequences': sample_sequences,
                'ref_sequences': ref_sequences,
                'outdir_alignments': outdir_alignments,
                'threads': threads,
                'gene_type': gene_type
            })
        else:
            logger.debug(f"{"[DEBUG]:":10} No reference sequences found for gene {gene_name}")
    
    if not alignment_tasks:
        logger.warning(f"{"[WARNING]:":10} No per-gene alignment tasks found")
        return
    
    logger.debug(f"{"[DEBUG]:":10} Found {len(alignment_tasks)} per-gene alignment tasks")
    
    # Process alignments with multiprocessing and progress bar
    logger.info(f"{"[INFO]:":10} Processing {len(alignment_tasks)} per-gene alignment tasks with {pool_size} processes and {threads} threads per process")
    with ProcessPoolExecutor(max_workers=pool_size) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(create_single_per_gene_alignment, task): task['gene_name']
            for task in alignment_tasks
        }
        
        # Process results with progress bar
        for future in tqdm(as_completed(future_to_task), total=len(future_to_task), desc=f"{"[INFO]:":10} {"Generating per-gene alignments":<40}", file=sys.stdout):
            
            task_id = future_to_task[future]

            success, result = future.result()
            
            if success:
                logger.debug(f"{"[DEBUG]:":10} {task_id}: {result}")
            else:
                utils.log_manager.handle_error(result[0], result[1], "generate_per_gene_alignments()", task_id)

    logger.debug(f"{"[DEBUG]:":10} Per-gene alignment generation complete")


def align_genes(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads, refs_order):
    """
    Generate alignments for annotated genes with reference genes.
    
    Args:
        all_sample_results (dict): Results from gene checking (contains sequence_gene_data for each sample)
        ref_gene_seqrecords (dict): Dictionary of reference gene sequences
        output_directory (str): Directory to write the output files
        pool_size (int): Number of processes to use for multiprocessing
        threads (int): Number of threads to use for each process
        refs_order (list): List of orders to use for reference genes
    """
    
    logger.info(f"{"[INFO]:":10} Generating gene alignments...")
    
    generate_cds_alignments(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads)
    generate_rRNA_alignments(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads)
    generate_tRNA_alignments(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads)
    generate_per_gene_alignments(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads)
    
    logger.info(f"{"[INFO]:":10} Gene alignment generation complete")
    utils.log_separator(logger)


def detect_multi_sequence_fasta(fasta_file, logger=None):
    """
    Detect if a FASTA file contains multiple sequences.
    
    Args:
        fasta_file (str): Path to the FASTA file
        logger: Logger object (optional)
        
    Returns:
        tuple: (is_multi_sequence, sequence_count)
    """
    try:
        sequence_count = 0
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sequence_count += 1
                if sequence_count > 1:
                    return True, sequence_count
        
        return False, sequence_count
        
    except Exception as e:
        if logger:
            logger.error(f"{"[ERROR]:":10} Failed to detect sequences in {fasta_file}: {str(e)}")
        return False, 0


def split_multi_sequence_fasta(fasta_file, output_dir, filename_prefix, logger=None):
    """
    Split a FASTA file containing multiple sequences into separate files.
    
    Args:
        fasta_file (str): Path to the input FASTA file
        output_dir (str): Directory to write individual sequence files
        filename_prefix (str): Base filename prefix for the sequences
        logger: Logger object (optional)
        
    Returns:
        list: List of paths to individual sequence files
    """
    sequence_files = []
    sequence_count = 0
    
    try:
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sequence_count += 1
                
                # Create a safe sequence name (replace spaces and special chars)
                safe_seq_name = str(record.id).replace(' ', '_').replace(':', '_').replace('|', '_')
                
                # Create individual sequence file
                seq_filename = f"{filename_prefix}_seq{sequence_count:03d}_{safe_seq_name}.fasta"
                seq_filepath = os.path.join(output_dir, seq_filename)
                
                # Write individual sequence to file
                with open(seq_filepath, 'w') as seq_handle:
                    SeqIO.write(record, seq_handle, 'fasta')
                
                sequence_files.append(seq_filepath)
                
        if logger:
            logger.debug(f"{"[DEBUG]:":10} Split {fasta_file} into {sequence_count} individual sequence files")
        return sequence_files
        
    except Exception as e:
        if logger:
            logger.error(f"{"[ERROR]:":10} Failed to split multi-sequence FASTA file {fasta_file}: {str(e)}")
        return []


def validate_linearisation_genes(linearise_genes, gene_synonyms, logger=None):
    """
    Validate that linearisation genes are present in the gene_synonyms.txt file.
    
    Args:
        linearise_genes (list): List of gene names to validate
        gene_synonyms (dict): Dictionary mapping gene names to standardized synonyms
        logger: Logger instance for logging messages
        
    Returns:
        list: List of validated gene names (RHS values from gene_synonyms.txt)
        
    Raises:
        ValueError: If any gene is not found in gene_synonyms.txt
    """
    validated_genes = []
    gene_synonyms_rhs_values = set(gene_synonyms.values())

    failed_genes = []
    
    for gene in linearise_genes:
        gene = gene.strip()
        if gene in gene_synonyms_rhs_values:
            validated_genes.append(gene)
        else:
            # Check if it's a LHS value that maps to a RHS value
            if gene in gene_synonyms:
                mapped_gene = gene_synonyms[gene]
                validated_genes.append(mapped_gene)
                if logger:
                    logger.debug(f"{"[DEBUG]:":10} Mapped linearisation gene: {gene} -> {mapped_gene}")
            else:
                failed_genes.append(gene)
             
    if failed_genes:
        logger.error(f"{"[ERROR]:":10} Failed to find the following linearisation genes in gene_synonyms.txt: {', '.join(failed_genes)}")
        utils.exit_program()
    
    logger.info(f"{"[INFO]:":10} Using linearisation genes: {', '.join(validated_genes)}")
    utils.log_separator(logger)
    
    return validated_genes


def linearise_genome_upstream_gene(gbk_file, fasta_file, output_dir, sample_name, linearise_genes=['psbA'], logger=None):
    """
    Linearise genome upstream of specified gene and write new fasta file.
    
    Args:
        gbk_file (str): Path to annotated GenBank file
        fasta_file (str): Path to original fasta file
        output_dir (str): Output directory for linearised fasta
        sample_name (str): Sample name for output file
        linearise_genes (list): List of gene names to try for linearisation (default: ['psbA'])
        logger: Logger instance for logging messages

    Returns:
        str: Path to the new linearised fasta file
    """
    try:
        # Read the annotated GenBank file
        records = utils.parse_genbank_file(gbk_file, logger)
        record = records[0]
        
        # Try each gene in the list until one is found
        gene_start = None
        found_gene = None
        
        for linearise_gene in linearise_genes:
            for feature in record.features:
                if feature.type == 'CDS' and 'gene' in feature.qualifiers:
                    if feature.qualifiers['gene'][0] == linearise_gene:
                        gene_start = feature.location.start
                        found_gene = linearise_gene
                        break
            if gene_start is not None:
                break
        
        if gene_start is None:
            gene_list_str = ', '.join(linearise_genes)
            logger.debug(f"{"[DEBUG]:":10} None of the linearisation genes [{gene_list_str}] found in {sample_name} record {record.id}, using original sequence instead")
            return fasta_file
        
        # Read original fasta sequence
        with open(fasta_file, 'r') as handle:
            fasta_record = next(SeqIO.parse(handle, 'fasta'))
        
        # Linearise upstream of the specified gene
        sequence = fasta_record.seq
        linearised_sequence = sequence[gene_start:] + sequence[:gene_start]
        
        # Create new fasta record
        linearised_record = SeqRecord(
            seq=linearised_sequence,
            id=fasta_record.id,
            description=f"Linearised upstream of {found_gene} (original position: {gene_start})"
        )
        
        # Write linearised fasta file
        linearised_fasta = os.path.join(output_dir, f"{sample_name}.round2.fasta")
        with open(linearised_fasta, 'w') as handle:
            SeqIO.write(linearised_record, handle, 'fasta')
        
        logger.debug(f"{"[INFO]:":10} Linearised {sample_name} upstream of {found_gene} (position {gene_start})")
        return linearised_fasta
        
    except Exception as e:
        logger.error(f"{"[ERROR]:":10} Failed to linearise {sample_name}: {str(e)}")
        return fasta_file


def process_single_sequence(fasta_file, output_dir, sequence_name, chloe_project_dir, linearise_genes, 
                            metadata_dict, original_fasta_name, logger=None):
    """
    Process a single sequence file with Chloë annotation.
    
    Args:
        fasta_file (str): Path to the sequence FASTA file
        output_dir (str): Output directory for this sequence
        sequence_name (str): Name for this sequence
        chloe_project_dir (str): Chloë project directory
        chloe_script_path (str): Chloë script path
        linearise_genes (list): List of genes to try for linearisation
        metadata_dict (dict): Metadata dictionary
        original_fasta_name (str): Original multi-sequence filename for metadata lookup
        logger: Logger instance for logging messages
        
    Returns:
        dict: Dictionary with GenBank and GFF file paths for this sequence
    """
    # Define output file paths for this sequence
    output_gbk = os.path.join(output_dir, f"{sequence_name}.chloe.gbk")
    output_gff = os.path.join(output_dir, f"{sequence_name}.chloe.gff")
    output_gbk_original = os.path.join(output_dir, f"{sequence_name}.round1.chloe.gbk")
    output_gff_original = os.path.join(output_dir, f"{sequence_name}.round1.chloe.gff")
    output_gbk_linearised = os.path.join(output_dir, f"{sequence_name}.round2.chloe.gbk")
    output_gff_linearised = os.path.join(output_dir, f"{sequence_name}.round2.chloe.gff")
    output_fasta_linearised = os.path.join(output_dir, f"{sequence_name}.round2.fasta")

    # Remove any existing non-renamed files
    if os.path.exists(output_gbk):
        os.remove(output_gbk)
    if os.path.exists(output_gff):
        os.remove(output_gff)

    # Determine what needs to be done
    has_original_annotation = (utils.file_exists_and_not_empty(output_gbk_original) and 
                               utils.file_exists_and_not_empty(output_gff_original))
    has_linearised_annotation_and_fasta = (utils.file_exists_and_not_empty(output_gbk_linearised) and 
                                           utils.file_exists_and_not_empty(output_gff_linearised) and 
                                           utils.file_exists_and_not_empty(output_fasta_linearised))
    
    if has_linearised_annotation_and_fasta:
        logger.debug(f"{"[INFO]:":10} {sequence_name} already completely processed with linearised sequence")
    elif has_original_annotation:
        logger.debug(f"{"[INFO]:":10} {sequence_name} needs linearisation and re-annotation")
    else:
        logger.debug(f"{"[INFO]:":10} {sequence_name} needs initial annotation")

    chloe_project = f'--project={chloe_project_dir}'
    chloe_script = os.path.join(chloe_project_dir, 'chloe.jl')
    
    # Run chloe annotate command
    cmd = [
        'julia',
        chloe_project,
        chloe_script,
        'annotate',
        '--no-filter',
        '--no-transform',
        '--gbk',
        '--output',
        output_dir,
        fasta_file
    ]

    # Execute annotation workflow
    if not has_original_annotation:
        # Do initial annotation
        logger.debug(f"{"[INFO]:":10} Performing initial annotation for {sequence_name}")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        # Rename original annotation files to preserve original annotation
        if os.path.exists(output_gbk):
            os.rename(output_gbk, output_gbk_original)
        if os.path.exists(output_gff):
            os.rename(output_gff, output_gff_original)

        has_original_annotation = True

    # Check metadata to determine if linearisation should be skipped
    should_linearise = True
    successful_linearisation = False
    if metadata_dict:
        # Use original multi-sequence filename for metadata lookup
        if original_fasta_name in metadata_dict:
            topology = metadata_dict[original_fasta_name].get('linear_or_circular', 'circular')
            if topology == 'linear':
                should_linearise = False
                logger.debug(f"{"[INFO]:":10} Skipping linearisation for {sequence_name} - already marked as linear in metadata")

    # Linearise and re-annotate (if we have original annotation but no linearised and sample should be linearised)
    if has_original_annotation and not has_linearised_annotation_and_fasta and should_linearise:
        gene_list_str = ', '.join(linearise_genes)
        logger.debug(f"{"[INFO]:":10} Linearising {sequence_name} upstream of genes: {gene_list_str}...")
        
        # Linearise the genome upstream of the specified gene
        output_fasta_linearised = linearise_genome_upstream_gene(
            output_gbk_original, 
            fasta_file, 
            output_dir, 
            sequence_name,
            linearise_genes,
            logger
        )

        if output_fasta_linearised != fasta_file:
            successful_linearisation = True
            # Re-run chloe with linearised fasta
            cmd[9] = output_fasta_linearised  # Replace the fasta file in the command

            if logger:
                logger.debug(f"{"[DEBUG]:":10} Re-annotating {sequence_name} with linearised sequence")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            if logger:
                logger.debug(f"{"[DEBUG]:":10} Re-annotation complete for {sequence_name}")

    elif has_original_annotation and has_linearised_annotation_and_fasta:
        successful_linearisation = True

    # Return the appropriate file paths
    if should_linearise and successful_linearisation:
        return {
            'gbk': output_gbk_linearised,
            'gff': output_gff_linearised,
            # 'fasta': output_fasta_linearised
        }
    else:
        return {
            'gbk': output_gbk_original,
            'gff': output_gff_original,
            # 'fasta': fasta_file
        }


def annotate_genomes(genome_fasta_dir, output_directory, chloe_project_dir,
                     linearise_genes=['psbA'],  metadata_dict=None, pool_size=1, log_queue=None):
    """
    Annotate genome fasta files using chloe annotate command.
    
    Args:
        genome_fasta_dir (str): Path to directory containing genome fasta files
        output_directory (str): Path to directory to write annotated genomes
        chloe_project_dir (str, optional): Path to the chloe project directory
        linearise_genes (list, optional): List of genes to try for linearisation (default: ['psbA'])
        metadata_dict (dict, optional): Metadata dictionary mapping fasta filenames to metadata.
                                      Used to check linear_or_circular status to skip linearisation
                                      for samples already marked as linear.
        pool_size (int, optional): Number of processes to use for multiprocessing. Defaults to 1.
        log_queue (queue.Queue, optional): Multiprocessing-safe queue for logging messages.
        
    Returns:
        dict: Dictionary mapping sample names to their annotation file paths
        
    Raises: 
        FileNotFoundError: If the directory doesn't exist
        Exception: For other errors during annotation
    """
    # Check if directory exists
    if not os.path.exists(genome_fasta_dir):
        logger.error(f"{"[ERROR]:":10} Directory not found: {genome_fasta_dir}")
        utils.exit_program()

    # Create a directory for the annotated genomes
    annotated_genomes_dir = os.path.join(output_directory, "01_annotated_genomes")
    os.makedirs(annotated_genomes_dir, exist_ok=True)

    # Find all fasta files in the input directory
    fasta_files = glob.glob(os.path.join(genome_fasta_dir, "*.fasta")) + \
                 glob.glob(os.path.join(genome_fasta_dir, "*.fa")) + \
                 glob.glob(os.path.join(genome_fasta_dir, "*.fas"))

    if not fasta_files:
        logger.error(f"{"[ERROR]:":10} No fasta files found in {genome_fasta_dir}")
        utils.exit_program()
    
    # Print fasta files to annotate
    logger.info(f"{"[INFO]:":10} Found {len(fasta_files)} fasta files to annotate:")
    for fasta_file in fasta_files:
        logger.info(f"{" ":10} {os.path.basename(fasta_file)}")
    utils.log_separator(logger)
    time.sleep(0.1)  # As tqdm output is not written to the queue, we need to sleep to ensure the logger messages have time to be written 

    annotated_genomes = dict({})
    
    # Process files with multiprocessing
    logger.info(f"{"[INFO]:":10} Processing {len(fasta_files)} FASTA files with {pool_size} processes")
    
    with ProcessPoolExecutor(max_workers=pool_size) as executor:
        # Submit all tasks
        future_to_file = {
            executor.submit(process_single_fasta_file, fasta_file, annotated_genomes_dir, chloe_project_dir,
                            linearise_genes, metadata_dict, log_queue): fasta_file
            for fasta_file in fasta_files
        }
        
        # Process results with progress bar
        for future in tqdm(as_completed(future_to_file), total=len(future_to_file),
                           desc=f"{"[INFO]:":10} {"Annotating genomes":<20}", file=sys.stdout):
            fasta_file = future_to_file[future]
            
            success, result = future.result()
                
            if success:
                sample_name, sample_data = result
                annotated_genomes[sample_name] = sample_data
                logger.debug(f"{"[DEBUG]:":10} Successfully processed {sample_name}")
            else:
                utils.log_manager.handle_error(result[0], result[1], 'annotate_genomes()', fasta_file)
                                               
    utils.log_separator(logger)
    time.sleep(0.1)

    return annotated_genomes


def process_single_fasta_file(fasta_file, annotated_genomes_dir, chloe_project_dir, linearise_genes,
                              metadata_dict, log_queue=None):
    """
    Process a single FASTA file with annotation (worker function for multiprocessing).
    
    Args:
        fasta_file (str): Path to the FASTA file to process
        annotated_genomes_dir (str): Directory for annotated genomes
        chloe_project_dir (str): Chloë project directory
        chloe_script_path (str): Chloë script path
        linearise_genes (list): List of genes to try for linearisation
        metadata_dict (dict): Metadata dictionary
        log_queue (queue.Queue, optional): Multiprocessing-safe queue for logging
        
    Returns:
        tuple: (success, result) where success is bool and result is either (sample_name, sample_data) or error info
    """
    try:
        worker_logger = utils.setup_worker_logger(__name__, log_queue)
        
        # Get the basename and filename prefix
        input_basename = os.path.basename(fasta_file)
        filename_prefix = os.path.splitext(input_basename)[0]

        # Remove additional extensions if present 
        filename_prefix = os.path.splitext(input_basename)[0]
        filename_prefix_no_dots = filename_prefix.replace('.', '_')

        # Check if this is a multi-sequence FASTA file
        is_multi_sequence, sequence_count = detect_multi_sequence_fasta(fasta_file, worker_logger)
        
        if is_multi_sequence:
            worker_logger.debug(f"{"[DEBUG]:":10} Multi-sequence FASTA detected: {input_basename} ({sequence_count} sequences)")
            
            # Create a subdirectory for this multi-sequence sample
            output_sample_dir = os.path.join(annotated_genomes_dir, filename_prefix_no_dots)
            os.makedirs(output_sample_dir, exist_ok=True)
            
            # Split the multi-sequence FASTA into individual files
            individual_sequence_files = split_multi_sequence_fasta(fasta_file, output_sample_dir, filename_prefix,
                                                                   worker_logger)
            
            if not individual_sequence_files:
                worker_logger.error(f"{"[ERROR]:":10} Failed to split multi-sequence FASTA file: {input_basename}")
                return False, (f'Failed to split multi-sequence FASTA file: {input_basename}', traceback.format_exc())
            
            # Initialize dict for multi-sequence sample data
            sample_data = {
                'input_filename': input_basename,
                'is_multi_sequence': True,
                'sequence_count': sequence_count,
                'sequences': {}
            }
            
            # Process each individual sequence
            for seq_idx, seq_file in enumerate(individual_sequence_files, 1):
                seq_basename = os.path.basename(seq_file)
                seq_name = os.path.splitext(seq_basename)[0]
                
                worker_logger.debug(f"{"[DEBUG]:":10} Processing sequence {seq_idx}/{sequence_count}: {seq_name}")
                
                # Process this individual sequence
                seq_result = process_single_sequence(
                    seq_file, 
                    output_sample_dir, 
                    seq_name, 
                    chloe_project_dir, 
                    linearise_genes, 
                    metadata_dict, 
                    input_basename,  # Original multi-sequence filename for metadata lookup
                    worker_logger
                )
                
                if seq_result:
                    sample_data['sequences'][seq_name] = seq_result
            
            # Store the multi-sequence sample info
            # For multi-sequence, collect all individual sequence files
            if sample_data['sequences']:
                # Store all individual sequence files
                all_gbk_files = []
                all_gff_files = []
                
                for seq_name, seq_result in sample_data['sequences'].items():
                    all_gbk_files.append(seq_result['gbk'])
                    all_gff_files.append(seq_result['gff'])
                
                # Store all files in the main dictionary
                sample_data['gbk'] = all_gbk_files
                sample_data['gff'] = all_gff_files
            
            return True, (filename_prefix_no_dots, sample_data)
        
        else:
            # Single sequence FASTA - process as before
            worker_logger.debug(f"{"[DEBUG]:":10} Single sequence FASTA: {input_basename}")
            
            # Initialize dict for single sequence sample data
            sample_data = {
                'input_filename': input_basename,
                'is_multi_sequence': False,
                'sequence_count': 1
            }
            
            # Process single sequence
            result = process_single_sequence(
                fasta_file, 
                os.path.join(annotated_genomes_dir, filename_prefix_no_dots), 
                filename_prefix, 
                chloe_project_dir, 
                linearise_genes, 
                metadata_dict, 
                input_basename,
                worker_logger
            )
            
            if result:
                sample_data.update(result)
            
            return True, (filename_prefix_no_dots, sample_data)
            
    except Exception as e:
        return False, (e, traceback.format_exc())


def convert_gbk_to_embl(annotated_genomes_dict, output_directory, metadata_dict=None):
    """
    Convert GenBank files to EMBL format with sample-specific metadata.
    
    Args:
        annotated_genomes_dict (dict): Dictionary of annotated genome information
        output_directory (str): Output directory for EMBL files
        metadata_dict (dict, optional): Dictionary mapping input_filename to metadata containing:
            - project_id: Project identifier
            - locus_tag: Locus tag prefix
            - genus_species: Genus and species name
            - linear_or_circular: Genome topology (linear or circular)
    
    Returns:
        None
    """

    logger.info(f"{"[INFO]:":10} Converting GenBank files to EMBL format...")

    # Create EMBL output directory
    embl_output_dir = os.path.join(output_directory, '02_embl_files')
    os.makedirs(embl_output_dir, exist_ok=True)

    #Create ENA submission directory
    ena_submission_dir = os.path.join(embl_output_dir, 'ena_submission_embl_files')
    os.makedirs(ena_submission_dir, exist_ok=True)
    
    # Track metadata usage
    samples_with_metadata = 0
    
    for sample_name, genome_info in annotated_genomes_dict.items():

        if 'error' in genome_info:
            logger.warning(f"{"[WARNING]:":10} Skipping {sample_name} - has error: {genome_info['error']}")
            continue
            
        gbk_files = genome_info.get('gbk')
        
        # Handle both single files and lists of files (for multi-sequence samples)
        if isinstance(gbk_files, list):
            # Multi-sequence sample
            if not gbk_files or not all(os.path.exists(f) for f in gbk_files):
                logger.warning(f"{"[WARNING]:":10} One or more GenBank files not found for {sample_name}")
                continue
        else:
            # Single sequence sample
            if not gbk_files or not os.path.exists(gbk_files):
                logger.warning(f"{"[WARNING]:":10} GenBank file not found for {sample_name}")
                continue
            # Convert to lists for consistent processing
            gbk_files = [gbk_files]
        
        # Get fasta filename for metadata lookup
        input_filename = genome_info['input_filename']
        
        # Get metadata for this sample (all samples must be in metadata by this point if running main() for subcommand `annotate_and_check`)
        if metadata_dict and input_filename in metadata_dict:
            sample_metadata = metadata_dict[input_filename]
            samples_with_metadata += 1
            logger.debug(f"{"[DEBUG]:":10} Using metadata for {sample_name} (fasta: {input_filename})")
 
            # Extract some metadata values here
            locus_tag_prefix = sample_metadata['locus_tag']
            genus_species = sample_metadata['genus_species']
        else:
            logger.debug(f"{"[DEBUG]:":10} No metadata found for {sample_name} (fasta: {input_filename})")
            locus_tag_prefix = 'DEFAULT_TAG'
            genus_species = 'Unknown species'
            sample_metadata = None
            
        # Process each file in the lists and collect ENA EMBL files for concatenation
        ena_embl_files = []
        
        for file_idx, gbk_file in enumerate(gbk_files):
            try:
                # For multi-sequence samples, create sequence-specific names
                if len(gbk_files) > 1:
                    seq_suffix = f"_seq{file_idx+1:03d}"
                    embl_sample_name = f"{sample_name}{seq_suffix}"
                else:
                    embl_sample_name = sample_name
                
                # Read GenBank file   
                records = utils.parse_genbank_file(gbk_file, logger)
                record = records[0]
                record_seq_length = len(record.seq)

                # Add source feature at the beginning of features list
                source_location = FeatureLocation(0, len(record.seq))
                source_feature = SeqFeature(
                    location=source_location,
                    type='source',
                    qualifiers={
                        'organism': [genus_species],
                        'mol_type': ['genomic DNA']
                    }
                )
                record.features.insert(0, source_feature)

                # Remove all existing locus_tags
                for feature in record.features:
                    if 'locus_tag' in feature.qualifiers:
                        del feature.qualifiers['locus_tag']
                
                # Initialize locus counter before any special-case tagging
                locus_counter = 1

                # Special handling for rps12:
                # - Keep the two rps12 CDS features
                # - Create a corresponding gene feature for each CDS with identical coordinates
                # - Remove all other rps12 features (e.g., gene or intron features that don't match the CDS)
                try:
                    rps12_cds = []
                    for f in record.features:
                        if f.type == 'CDS' and 'gene' in f.qualifiers and f.qualifiers['gene'][0] == 'rps12':
                            rps12_cds.append(f)

                    if rps12_cds:
                        # Remove all non-CDS rps12 features
                        new_features = []
                        for f in record.features:
                            gene_name = None
                            if 'gene' in f.qualifiers:
                                gene_name = f.qualifiers['gene'][0]
                            elif 'product' in f.qualifiers:
                                gene_name = f.qualifiers['product'][0]
                            if gene_name == 'rps12' and f.type != 'CDS':
                                continue
                            new_features.append(f)

                        # Insert a gene feature immediately before each rps12 CDS
                        for cds in rps12_cds:
                            try:
                                insert_at = new_features.index(cds)
                            except ValueError:
                                insert_at = len(new_features)
                            rps12_gene_feature = SeqFeature(
                                location=cds.location,
                                type='gene',
                                qualifiers={'gene': ['rps12']}
                            )
                            # Assign a shared locus_tag for the CDS and its paired gene feature
                            pair_tag = f"{locus_tag_prefix}_LOCUS{locus_counter}"
                            locus_counter += 1
                            rps12_gene_feature.qualifiers['locus_tag'] = [pair_tag]
                            cds.qualifiers['locus_tag'] = [pair_tag]
                            new_features.insert(insert_at, rps12_gene_feature)

                        record.features = new_features
                except Exception:
                    # Fail-safe: proceed without special handling if any unexpected structure is encountered
                    pass

                # Modify features for EMBL output
                
                # First pass: identify gene groups and assign locus_tags in order
                gene_occurrence_map = {}  # gene_name -> occurrence_count
                gene_group_locus_map = {}  # (gene_name, occurrence) -> locus_tag
                feature_to_gene_group = {}  # feature_index -> (gene_name, occurrence)
                
                for idx, feature in enumerate(record.features):
                    if feature.type in ['gene', 'CDS', 'tRNA', 'rRNA', 'intron']:
                        # Skip if already has a locus_tag (e.g., from rps12 special handling)
                        if 'locus_tag' in feature.qualifiers:
                            continue

                        # Get gene name
                        gene_name = None
                        if 'gene' in feature.qualifiers:
                            gene_name = feature.qualifiers['gene'][0]
      
                        if gene_name:
                            # Determine which occurrence of this gene we're dealing with
                            if feature.type == 'gene':
                                # New gene occurrence
                                gene_occurrence_map[gene_name] = gene_occurrence_map.get(gene_name, 0) + 1
                                current_occurrence = gene_occurrence_map[gene_name]
                                
                                # Assign new locus tag for this gene occurrence
                                gene_group = (gene_name, current_occurrence)
                                gene_group_locus_map[gene_group] = f"{locus_tag_prefix}_LOCUS{locus_counter}"
                                locus_counter += 1
                            else:
                                # Feature belongs to the most recent occurrence of this gene
                                current_occurrence = gene_occurrence_map.get(gene_name, 1)
                                gene_group = (gene_name, current_occurrence)
                                
                                # Create locus tag if this is the first feature for a gene without a gene feature
                                if gene_group not in gene_group_locus_map:
                                    gene_occurrence_map[gene_name] = current_occurrence
                                    gene_group_locus_map[gene_group] = f"{locus_tag_prefix}_LOCUS{locus_counter}"
                                    locus_counter += 1
                            
                            feature_to_gene_group[idx] = gene_group
                        
                    elif feature.type == 'repeat_region':
                        feature.qualifiers['locus_tag'] = [f"{locus_tag_prefix}_LOCUS{locus_counter}"]
                        locus_counter += 1

                # Second pass: apply locus tags and other modifications
                for idx, feature in enumerate(record.features):
                    # Remove /ID= and /parent= and /name qualifiers if present (not used in EMBL format)
                    if 'ID' in feature.qualifiers:
                        del feature.qualifiers['ID']
                    if 'parent' in feature.qualifiers:
                        del feature.qualifiers['parent']
                    if 'name' in feature.qualifiers:
                        del feature.qualifiers['name']
                    
                    if feature.type in ['gene', 'CDS', 'tRNA', 'rRNA', 'intron']:
                        # Skip if already has a locus_tag (e.g., from rps12 special handling)
                        if 'locus_tag' not in feature.qualifiers:
                            if idx in feature_to_gene_group:
                                gene_group = feature_to_gene_group[idx]
                                if gene_group in gene_group_locus_map:
                                    feature.qualifiers['locus_tag'] = [gene_group_locus_map[gene_group]]
                    
                    # Additional modifications for CDS features
                    if feature.type == 'CDS':
                        # Add /transl_table=11 to all CDS features
                        feature.qualifiers['transl_table'] = ['11']

                        # Add /codon_start=1 to all CDS features
                        feature.qualifiers['codon_start'] = ['1']
                        
                        # Add /trans_splicing to rps12 genes
                        if 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] == 'rps12':
                            feature.qualifiers['trans_splicing'] = ['']
                
                # Convert to EMBL format
                embl_filename = f"{embl_sample_name}.embl"
                embl_filepath = os.path.join(embl_output_dir, embl_filename)
                
                # Write EMBL file first pass
                with open(embl_filepath, 'w') as handle:
                    SeqIO.write(record, handle, 'embl')

                # Convert EMBL file to ENA template format
                ena_embl_filepath = convert_embl_to_ena_template(embl_filepath, sample_metadata, record_seq_length, record.id, ena_submission_dir)
                
                # Store the ENA EMBL file path for later concatenation
                ena_embl_files.append(ena_embl_filepath)
                
                logger.debug(f"{"[DEBUG]:":10} Converted {embl_sample_name} to EMBL format: {ena_embl_filepath}")
                
            except Exception as e:
                logger.error(f"{"[ERROR]:":10} Failed to convert {embl_sample_name if 'embl_sample_name' in locals() else sample_name} to EMBL format: {str(e)}\n{traceback.format_exc()}")
                continue
        
        # Concatenate all ENA EMBL files for this sample into a single file
        if ena_embl_files:
            try:
                # Create the concatenated file name
                if len(ena_embl_files) == 1:
                    # Single file - just rename it
                    final_ena_embl_filepath = ena_embl_files[0]
                else:
                    # Multiple files - concatenate them
                    final_ena_embl_filepath = os.path.join(ena_submission_dir, f"{sample_name}.ena.embl")
                    
                    with open(final_ena_embl_filepath, 'w') as outfile:
                        for i, ena_embl_file in enumerate(ena_embl_files):
                            # Read and write the content of each file
                            with open(ena_embl_file, 'r') as infile:
                                content = infile.read()
                                outfile.write(content)
                            
                            # Add separator '//' between records (except after the last one)
                            if i < len(ena_embl_files) - 1:
                                outfile.write('//\n')
                    
                    # Remove individual ENA EMBL files after concatenation
                    for ena_embl_file in ena_embl_files:
                        try:
                            os.remove(ena_embl_file)
                            logger.debug(f"{"[DEBUG]:":10} Removed individual ENA EMBL file: {ena_embl_file} from {ena_submission_dir}")
                        except Exception as e:
                            logger.warning(f"{"[WARNING]:":10} Failed to remove individual ENA EMBL file {ena_embl_file}: {e}")
                    
                    logger.info(f"{"[INFO]:":10} Concatenated {len(ena_embl_files)} ENA EMBL files for {sample_name} into: {final_ena_embl_filepath}")
                
            except Exception as e:
                logger.error(f"{"[ERROR]:":10} Failed to concatenate ENA EMBL files for {sample_name}: {str(e)}")
                continue
    
    logger.info(f"{"[INFO]:":10} EMBL conversion complete. Files written to: {embl_output_dir}")
    
    # Log metadata usage summary
    logger.debug(f"{"[DEBUG]:":10} Successfully processed {samples_with_metadata} samples with metadata")
    
    utils.log_separator(logger)


def convert_embl_to_ena_template(embl_filepath, sample_metadata, record_seq_length, record_id, ena_submission_dir):
    """
    Convert EMBL file to ENA template format.
    
    This function takes an EMBL file and converts it to the European Nucleotide Archive (ENA)
    template format, incorporating sample metadata for submission.
    
    Args:
        embl_filepath (str): Path to the input EMBL file
        sample_metadata (dict): Dictionary containing sample metadata with keys:
            - locus_tag (str): Locus tag prefix for the genome
            - project_id (str): Project identifier
            - genus_species (str): Genus and species name
            - linear_or_circular (str): Whether genome is linear or circular
        record_seq_length (int): Length of the sequence in the record
        record_id (str): ID of the record
    
    Returns:
        None: Modifies the EMBL file in place or creates a new ENA template file
        
    Raises:
        FileNotFoundError: If the EMBL file does not exist
        ValueError: If required metadata fields are missing
        Exception: For other processing errors
    """

    logger.debug(f"{"[DEBUG]:":10} Converting EMBL file to ENA template: {embl_filepath}")
    logger.debug(f"{"[DEBUG]:":10} Using metadata: {sample_metadata}")

    lines_to_write = []
    space_line_template = 'XX'

    id_line_template = f"ID   XXX; XXX; {sample_metadata['linear_or_circular'] if sample_metadata else 'unknown'}; XXX; XXX; XXX; XXX."
    lines_to_write.append(id_line_template)
    lines_to_write.append(space_line_template)

    ac_line_template = f"AC   XXX;"
    lines_to_write.append(ac_line_template)
    lines_to_write.append(space_line_template)

    ac_line_template_2 = f"AC * _{record_id};"
    lines_to_write.append(ac_line_template_2)
    lines_to_write.append(space_line_template)

    pr_line_template = f"PR   Project:{sample_metadata['project_id'] if sample_metadata else 'UNKNOWN_PROJECT'};"
    lines_to_write.append(pr_line_template)
    lines_to_write.append(space_line_template)

    de_line_template = f"DE   XXX"
    lines_to_write.append(de_line_template)
    lines_to_write.append(space_line_template)

    rn_line_template = f"RN   [1]"
    lines_to_write.append(rn_line_template)

    rp_line_template = f"RP   1-{record_seq_length}"
    lines_to_write.append(rp_line_template)

    ra_line_template = f"RA   XXX;"
    lines_to_write.append(ra_line_template)

    rt_line_template = f"RT   ;"
    lines_to_write.append(rt_line_template)
    
    rl_line_template = f"RL   Submitted (DD-MMM-YYYY) to the INSDC."
    lines_to_write.append(rl_line_template)
    lines_to_write.append(space_line_template)

    # Get FH and FT lines

    fh_lines = []
    ft_lines = []
    sq_lines = []
    nucleotide_lines = []

    with open(embl_filepath, 'r') as handle:
        for line in handle:
            if line.startswith('FH'):
                fh_lines.append(line.strip())
            elif line.startswith('FT'):
                ft_lines.append(line.strip())
            elif line.startswith('SQ'):
                sq_lines.append(line.strip())
            elif line.startswith('     '):
                nucleotide_lines.append(line.rstrip())

    lines_to_write.extend(fh_lines)
    lines_to_write.extend(ft_lines)
    lines_to_write.append(space_line_template)
    lines_to_write.extend(sq_lines)
    lines_to_write.extend(nucleotide_lines)

    # Write ENA template file
    base_name = os.path.basename(embl_filepath)
    ena_template_filepath = os.path.join(ena_submission_dir, base_name.replace('.embl', '.ena.embl'))
    with open(ena_template_filepath, 'w') as handle:
        for line in lines_to_write:
            handle.write(line + '\n')   

    logger.debug(f"{"[DEBUG]:":10} ENA template file written to: {ena_template_filepath}")
    
    return ena_template_filepath
   
    
def detect_blast_database(directory_path, logger=None):
    """
    Detect BLAST database files in a directory and return the database path.
    
    Args:
        directory_path (str): Path to directory containing BLAST database files
        logger: Logger instance for logging messages
        
    Returns:
        str: Full path to the BLAST database (without extension) if found, None otherwise
    """
    if not os.path.exists(directory_path):
        if logger:
            logger.error(f"{"[ERROR]:":10} Directory does not exist: {directory_path}")
        return None
    
    if not os.path.isdir(directory_path):
        if logger:
            logger.error(f"{"[ERROR]:":10} Path is not a directory: {directory_path}")
        return None
    
    # Look for BLAST database files (common extensions: .nhr, .nin, .nsq for nucleotide)
    blast_extensions = ['.nhr', '.nin', '.nsq', '.ndb', '.not', '.ntf', '.nto']
    database_names = set()
    
    for filename in os.listdir(directory_path):
        for ext in blast_extensions:
            if filename.endswith(ext):
                # Extract database name by removing extension
                db_name = filename[:-len(ext)]
                database_names.add(db_name)
                break
    
    if not database_names:
        if logger:
            logger.error(f"{"[ERROR]:":10} No BLAST database files found in directory: {directory_path}")
            logger.error(f"{"":10} Expected files with extensions: {', '.join(blast_extensions)}")
        return None
    
    if len(database_names) > 1:
        if logger:
            logger.warning(f"{"[WARNING]:":10} Multiple BLAST databases found in directory: {directory_path}")
            logger.warning(f"{"":10} Databases: {', '.join(sorted(database_names))}")
            logger.warning(f"{"":10} Using the first one: {sorted(database_names)[0]}")
    
    # Use the first (or only) database found
    db_name = sorted(database_names)[0]
    db_path = os.path.join(directory_path, db_name)
    
    if logger:
        logger.info(f"{"[INFO]:":10} Detected BLAST database: {db_name} in {directory_path}")
    
    return db_path


def query_intergenic_regions(annotated_genomes_dict, output_directory, min_intergenic_length=50, blast_evalue=1e-10, 
                             debug_intergenic=False, max_blast_hits=1, pool_size=1, threads=1, log_queue=None, custom_blast_db=None):
    """
    Query intergenic regions from annotated genomes.
    
    Args:
        annotated_genomes_dict (dict): Dictionary of annotated genome information
        output_directory (str): Output directory for results
        min_intergenic_length (int): Minimum length of intergenic region to analyze
        blast_evalue (float): BLAST E-value threshold
        debug_intergenic (bool): Whether to write debug FASTA files
        max_blast_hits (int): Maximum number of BLAST hits to retain per region
        pool_size (int): Number of processes to use for multiprocessing
        threads (int): Number of threads to use for each process
        log_queue (queue.Queue, optional): Multiprocessing-safe queue for logging
        custom_blast_db (str, optional): Custom BLAST database path. If None, uses default order_genomes_blastdb
    
    Returns:
        None
    """
    logger.info(f"{"[INFO]:":10} Querying intergenic regions from annotated genomes...")

    # Create output directory for intergenic analysis
    intergenic_output_dir = os.path.join(output_directory, "04_intergenic_analysis")
    os.makedirs(intergenic_output_dir, exist_ok=True)

    # Path to BLAST database
    if custom_blast_db:
        # Check if the provided path is a directory (containing BLAST database files)
        if os.path.isdir(custom_blast_db):
            # Detect BLAST database in the directory
            blast_db_path = detect_blast_database(custom_blast_db, logger)
            if blast_db_path is None:
                logger.error(f"{"[ERROR]:":10} Failed to detect BLAST database in directory: {custom_blast_db}")
                return
        else:
            # Assume it's a direct path to the BLAST database
            blast_db_path = custom_blast_db
            logger.info(f"{"[INFO]:":10} Using custom BLAST database: {blast_db_path}")
    else:
        # Use default BLAST database
        blast_db_path = os.path.join("plastid_annotation_validator", "data", "order_genomes_blastdb", "order_genomes_blastdb")
        
        # Try alternative paths if the first one doesn't work
        if not os.path.exists(f"{blast_db_path}.nhr"):
            # Try relative to current working directory
            blast_db_path = os.path.join("plastid_annotation_validator", "data", "order_genomes_blastdb", "order_genomes_blastdb")
            if not os.path.exists(f"{blast_db_path}.nhr"):
                # Try absolute path from package
                import plastid_annotation_validator
                package_dir = os.path.dirname(plastid_annotation_validator.__file__)
                blast_db_path = os.path.join(package_dir, "data", "order_genomes_blastdb", "order_genomes_blastdb")
    
    # Check if BLAST database exists
    if not os.path.exists(f"{blast_db_path}.nhr"):
        logger.error(f"{"[ERROR]:":10} BLAST database not found at {blast_db_path}")
        logger.error(f"{"[ERROR]:":10} Please ensure the BLAST database has been created")
        return

    # Prepare data for multiprocessing and check for existing reports
    sample_data_list = []
    skipped_samples = []
    existing_results = []
    
    for sample_name, genome_info in annotated_genomes_dict.items():
        if 'error' in genome_info:
            logger.warning(f"{"[WARNING]:":10} Skipping {sample_name} - has error: {genome_info['error']}")
            continue
        
        # Check if report already exists
        genome_report_file = os.path.join(intergenic_output_dir, f"{sample_name}_intergenic_blast_results.tsv")
        if os.path.exists(genome_report_file) and os.path.getsize(genome_report_file) > 0:
            logger.debug(f"{"[INFO]:":10} Skipping {sample_name} - intergenic report already exists: {genome_report_file}")
            skipped_samples.append(sample_name)
            
            # Load existing results for combined report
            try:
                existing_blast_results = load_existing_blast_results(genome_report_file, sample_name)
                existing_results.extend(existing_blast_results)
                logger.debug(f"{"[DEBUG]:":10} Loaded {len(existing_blast_results)} existing BLAST results for {sample_name}")
            except Exception as e:
                logger.warning(f"{"[WARNING]:":10} Failed to load existing results for {sample_name}: {str(e)}")
            continue
        
        # Handle both single and multi-sequence samples
        is_multi_sequence = genome_info.get('is_multi_sequence', False)
        
        if is_multi_sequence:
            # Multi-sequence sample - get lists of files
            gbk_files = genome_info.get('gbk', [])
            
            if not gbk_files:
                logger.warning(f"{"[WARNING]:":10} Skipping {sample_name} - missing GenBankfiles")
                continue
                
            # Create sample data for multi-sequence processing
            sample_data = {
                'sample_name': sample_name,
                'genome_info': genome_info,
                'gbk_files': gbk_files,
                'is_multi_sequence': True,
                'blast_db_path': blast_db_path,
                'min_intergenic_length': min_intergenic_length,
                'blast_evalue': blast_evalue,
                'debug_intergenic': debug_intergenic,
                'max_blast_hits': max_blast_hits,
                'threads': threads,
                'intergenic_output_dir': intergenic_output_dir
            }
        else:
            # Single sequence sample - get single file paths
            gbk_file = genome_info.get('gbk')
            
            if not gbk_file:
                logger.warning(f"{"[WARNING]:":10} Skipping {sample_name} - missing GenBank files")
                continue
                
            # Create sample data for single sequence processing
            sample_data = {
                'sample_name': sample_name,
                'genome_info': genome_info,
                'gbk_files': [gbk_file],
                'is_multi_sequence': False,
                'blast_db_path': blast_db_path,
                'min_intergenic_length': min_intergenic_length,
                'blast_evalue': blast_evalue,
                'debug_intergenic': debug_intergenic,
                'max_blast_hits': max_blast_hits,
                'threads': threads,
                'intergenic_output_dir': intergenic_output_dir
            }
        
        sample_data_list.append(sample_data)
    
    if not sample_data_list:
        if skipped_samples:
            logger.info(f"{"[INFO]:":10} All {len(skipped_samples)} samples already have intergenic reports, skipping processing")
            # Write combined report with existing results only
            if existing_results:
                combined_report_file = os.path.join(intergenic_output_dir, "combined_intergenic_blast_results.tsv")
                write_blast_report(existing_results, combined_report_file, "ALL_GENOMES", logger)
                logger.info(f"{"[INFO]:":10} Combined report written to: {combined_report_file}")
            return
        else:
            logger.warning(f"{"[WARNING]:":10} No valid samples found for intergenic analysis")
            return
    
    # Log processing summary
    if skipped_samples:
        logger.info(f"{"[INFO]:":10} Skipping {len(skipped_samples)} samples with existing reports: {', '.join(skipped_samples[:5])}{'...' if len(skipped_samples) > 5 else ''}")
    
    # Combined results for all genomes (start with existing results)
    all_results = existing_results.copy()
    
    # Process samples with multiprocessing and progress bar
    logger.info(f"{"[INFO]:":10} Processing {len(sample_data_list)} samples with {pool_size} processes")
    with ProcessPoolExecutor(max_workers=pool_size) as executor:
        # Submit all tasks
        future_to_sample = {
            executor.submit(process_single_genome_intergenic, sample_data, log_queue): sample_data['sample_name']
            for sample_data in sample_data_list
        }
        
        # Process results with progress bar
        for future in tqdm(as_completed(future_to_sample), total=len(future_to_sample), desc=f"{"[INFO]:":10} {"Processing intergenic regions":<40}", file=sys.stdout):
            sample_name = future_to_sample[future]
            
            try:
                success, result = future.result()
                
                if success:
                    blast_results, genome_report_file = result
                    all_results.extend(blast_results)
                    logger.debug(f"{"[DEBUG]:":10} Completed intergenic analysis for {sample_name}: {len(blast_results)} BLAST hits")
                else:
                    utils.log_manager.handle_error(result[0], "Intergenic region analysis", sample_name)
                    
            except Exception as e:
                logger.error(f"{"[ERROR]:":10} Error processing {sample_name}: {str(e)}")
                logger.error(traceback.format_exc())
                continue

    # Write combined report
    if all_results:
        combined_report_file = os.path.join(intergenic_output_dir, "combined_intergenic_blast_results.tsv")
        write_blast_report(all_results, combined_report_file, "ALL_GENOMES", logger)
        logger.info(f"{"[INFO]:":10} Combined report written to: {combined_report_file}")
    else:
        logger.warning(f"{"[WARNING]:":10} No BLAST results found for any genome")

    logger.info(f"{"[INFO]:":10} Intergenic region analysis complete. Results in: {intergenic_output_dir}")
    utils.log_separator(logger)


def extract_intergenic_regions(gbk_file, min_length=50, debug_intergenic=False, logger=None):
    """
    Extract intergenic regions from a GenBank file.
    
    Args:
        gbk_file (str): Path to GenBank file
        min_length (int): Minimum length of intergenic region to extract
    
    Returns:
        list: List of dictionaries containing intergenic region information
    """
    intergenic_regions = []
    
    # Parse GenBank file
    records = utils.parse_genbank_file(gbk_file, logger)
    record = records[0]
    
    # Get features with locations, but exclude features that are parts of genes or structural regions
    # We want intergenic regions between genes and within genes (introns)
    excluded_types = ['exon', 'source', 'LSC', 'SSC', 'inverted_repeat', 'intron', 'repeat_region']
    kept_features = [f for f in record.features if hasattr(f, 'location') and f.location is not None and f.type not in excluded_types]
    feature_types = set([f.type for f in kept_features])

    # Simple approach: Get coordinates for CDS regions (or tRNA coordinates for tRNAs)
    # Group features by gene name and get the appropriate coordinates
    gene_coordinates = []
    
    # Process each feature individually to handle multiple copies of the same gene
    for feature in kept_features:
        gene_name = get_gene_name(feature)
        
        # Handle CDS features
        if feature.type == 'CDS':
            if hasattr(feature.location, 'parts') and len(feature.location.parts) > 1:
                # Split gene - add each part
                for part in feature.location.parts:
                    gene_coordinates.append({
                        'gene_name': gene_name,
                        'start': part.start,
                        'end': part.end,
                        'type': 'CDS'
                    })
            else:
                # Single location
                gene_coordinates.append({
                    'gene_name': gene_name,
                    'start': feature.location.start,
                    'end': feature.location.end,
                    'type': 'CDS'
                })
        
        # Handle tRNA features (only if no CDS for this gene)
        elif feature.type == 'tRNA':
            # Check if we already have CDS coordinates for this gene
            existing_cds = [coord for coord in gene_coordinates if coord['gene_name'] == gene_name and coord['type'] == 'CDS']
            if not existing_cds:
                if hasattr(feature.location, 'parts') and len(feature.location.parts) > 1:
                    # Split gene - add each part
                    for part in feature.location.parts:
                        gene_coordinates.append({
                            'gene_name': gene_name,
                            'start': part.start,
                            'end': part.end,
                            'type': 'tRNA'
                        })
                else:
                    # Single location
                    gene_coordinates.append({
                        'gene_name': gene_name,
                        'start': feature.location.start,
                        'end': feature.location.end,
                        'type': 'tRNA'
                    })
        
        # Handle rRNA features (only if no CDS or tRNA for this gene)
        elif feature.type == 'rRNA':
            # Check if we already have CDS or tRNA coordinates for this gene
            existing_coords = [coord for coord in gene_coordinates if coord['gene_name'] == gene_name and coord['type'] in ['CDS', 'tRNA']]
            if not existing_coords:
                if hasattr(feature.location, 'parts') and len(feature.location.parts) > 1:
                    # Split gene - add each part
                    for part in feature.location.parts:
                        gene_coordinates.append({
                            'gene_name': gene_name,
                            'start': part.start,
                            'end': part.end,
                            'type': 'rRNA'
                        })
                else:
                    # Single location
                    gene_coordinates.append({
                        'gene_name': gene_name,
                        'start': feature.location.start,
                        'end': feature.location.end,
                        'type': 'rRNA'
                    })
        
        # Handle other gene types (misc_RNA, ncRNA, etc.) as fallback
        elif feature.type in ['misc_RNA', 'ncRNA']:
            # Check if we already have coordinates for this gene
            existing_coords = [coord for coord in gene_coordinates if coord['gene_name'] == gene_name]
            if not existing_coords:
                if hasattr(feature.location, 'parts') and len(feature.location.parts) > 1:
                    # Split gene - add each part
                    for part in feature.location.parts:
                        gene_coordinates.append({
                            'gene_name': gene_name,
                            'start': part.start,
                            'end': part.end,
                            'type': feature.type
                        })
                else:
                    # Single location
                    gene_coordinates.append({
                        'gene_name': gene_name,
                        'start': feature.location.start,
                        'end': feature.location.end,
                        'type': feature.type
                    })

    # Sort coordinates by start position
    gene_coordinates.sort(key=lambda x: x['start'])
    
    if logger:
        logger.debug(f"{"[DEBUG]:":10} Found {len(gene_coordinates)} gene coordinates for intergenic calculation")
    
    # Find intergenic regions between gene coordinates
    for i in range(len(gene_coordinates) - 1):
        current_coord = gene_coordinates[i]
        next_coord = gene_coordinates[i + 1]
        
        # Get end of current gene and start of next gene
        current_end = current_coord['end']
        next_start = next_coord['start']
        
        # Check if there's a gap between genes
        if next_start > current_end:
            intergenic_length = next_start - current_end
            
            if intergenic_length >= min_length:
                # Extract intergenic sequence
                intergenic_seq = str(record.seq[current_end:next_start])
                
                # Create intergenic region info
                intergenic_info = {
                    'start': current_end,
                    'end': next_start,
                    'length': intergenic_length,
                    'sequence': intergenic_seq,
                    'upstream_gene': current_coord['gene_name'],
                    'downstream_gene': next_coord['gene_name'],
                    'region_id': f"IG_{current_end}_{next_start}"
                }
                
                intergenic_regions.append(intergenic_info)
    
    # Check for intergenic region at the beginning (before first gene)
    if gene_coordinates:
        first_coord = gene_coordinates[0]
        if first_coord['start'] > 0:
            intergenic_length = first_coord['start']
            if intergenic_length >= min_length:
                intergenic_seq = str(record.seq[0:first_coord['start']])
                intergenic_info = {
                    'start': 0,
                    'end': first_coord['start'],
                    'length': intergenic_length,
                    'sequence': intergenic_seq,
                    'upstream_gene': 'START',
                    'downstream_gene': first_coord['gene_name'],
                    'region_id': f"IG_0_{first_coord['start']}"
                }
                intergenic_regions.append(intergenic_info)
    
    # Check for intergenic region at the end (after last gene)
    if gene_coordinates:
        last_coord = gene_coordinates[-1]
        if last_coord['end'] < len(record.seq):
            intergenic_length = len(record.seq) - last_coord['end']
            if intergenic_length >= min_length:
                intergenic_seq = str(record.seq[last_coord['end']:])
                intergenic_info = {
                    'start': last_coord['end'],
                    'end': len(record.seq),
                    'length': intergenic_length,
                    'sequence': intergenic_seq,
                    'upstream_gene': last_coord['gene_name'],
                    'downstream_gene': 'END',
                    'region_id': f"IG_{last_coord['end']}_{len(record.seq)}"
                }
                intergenic_regions.append(intergenic_info)
    
    # Write debug FASTA file if requested (for end regions)
    if debug_intergenic and intergenic_regions:
        debug_fasta_file = gbk_file.replace('.gbk', '_intergenic_debug.fasta')
        with open(debug_fasta_file, 'w') as f:
            for region in intergenic_regions:
                f.write(f">{region['region_id']}\n{region['sequence']}\n")
        if logger:
            logger.debug(f"{"[DEBUG]:":10} Wrote {len(intergenic_regions)} intergenic regions to debug FASTA: {debug_fasta_file}")
                    
    return intergenic_regions


def process_single_genome_intergenic(sample_data, log_queue=None):
    """
    Process intergenic regions for a single genome (worker function for multiprocessing).
    
    Args:
        sample_data (dict): Dictionary containing sample information and parameters
        
    Returns:
        tuple: (success, result) where success is bool and result is either (blast_results, genome_report_file) or error message
    """
    try:
        worker_logger = utils.setup_worker_logger(__name__, log_queue)
        
        sample_name = sample_data['sample_name']
        gbk_files = sample_data['gbk_files']
        is_multi_sequence = sample_data['is_multi_sequence']
        blast_db_path = sample_data['blast_db_path']
        min_intergenic_length = sample_data['min_intergenic_length']
        blast_evalue = sample_data['blast_evalue']
        debug_intergenic = sample_data['debug_intergenic']
        max_blast_hits = sample_data['max_blast_hits']
        threads = sample_data['threads']
        intergenic_output_dir = sample_data['intergenic_output_dir']
        
        # Validate files exist
        for gbk_file in gbk_files:
            if not os.path.exists(gbk_file):
                return False, f"GenBank file not found: {gbk_file}"
        
        # Process all sequences for this sample
        all_blast_results = []
        
        for i, gbk_file in enumerate(gbk_files):
            # Extract sequence name for multi-sequence samples
            if is_multi_sequence:
                seq_basename = os.path.basename(gbk_file)
                # Remove file extensions
                seq_name = seq_basename.replace('.chloe.gbk', '').replace('.round1.chloe.gbk', '').replace('.round2.chloe.gbk', '')
                # Extract just the sequence part using the known file naming pattern
                # The filename format from split_multi_sequence_fasta is: sample1_seq001_sequence1.chloe.gbk
                # We want to extract: seq001_sequence1
                # Use regex to match the pattern: sample_name + "_seq" + 3 digits + "_" + sequence_name
                pattern = rf"^{re.escape(sample_name)}_seq\d{{3}}_(.+)$"
                match = re.match(pattern, seq_name)
                if match:
                    # Extract the sequence part (everything after sample_name_seq001_)
                    seq_name = match.group(1)
                sequence_id = f"{sample_name}_{seq_name}"
            else:
                sequence_id = sample_name
            
            worker_logger.debug(f"{"[DEBUG]:":10} Processing sequence {i+1}/{len(gbk_files)}: {sequence_id}")

            # Extract intergenic regions
            intergenic_regions = extract_intergenic_regions(gbk_file, min_intergenic_length,
                                                            debug_intergenic, worker_logger)
            
            if not intergenic_regions:
                worker_logger.warning(f"{"[WARNING]:":10} No intergenic regions found for {sequence_id}")
                continue

            # BLAST intergenic regions - use sequence_id to include sequence information
            blast_results = blast_intergenic_regions(intergenic_regions, blast_db_path, blast_evalue, sequence_id,
                                                     max_blast_hits, threads, worker_logger)
            all_blast_results.extend(blast_results)
        
        if not all_blast_results:
            return False, (ValueError(f"No intergenic regions found for {sample_name}"), traceback.format_exc())

        # Write individual genome report (single report per sample)
        genome_report_file = os.path.join(intergenic_output_dir, f"{sample_name}_intergenic_blast_results.tsv")
        write_blast_report(all_blast_results, genome_report_file, sample_name, worker_logger)
        
        return True, (all_blast_results, genome_report_file)
        
    except Exception as e:
        return False, (e, traceback.format_exc())


def get_gene_name(feature):
    """
    Extract gene name from a feature.
    
    Args:
        feature: BioPython feature object
    
    Returns:
        str: Gene name or feature type if no gene name found
    """
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    elif 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    elif 'product' in feature.qualifiers:
        return feature.qualifiers['product'][0]
    else:
        return feature.type


def blast_intergenic_regions(intergenic_regions, blast_db_path, evalue_threshold, sample_name, max_blast_hits=1, threads=1, logger=None):
    """
    BLAST intergenic regions against the database.
    
    Args:
        intergenic_regions (list): List of intergenic region dictionaries
        blast_db_path (str): Path to BLAST database
        evalue_threshold (float): E-value threshold
        sample_name (str): Name of the sample
    
    Returns:
        list: List of BLAST results
    """
    blast_results = []
    
    for region in intergenic_regions:
        
        # Create temporary FASTA file for the intergenic region
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
            temp_fasta.write(f">{region['region_id']}\n{region['sequence']}\n")
            temp_fasta_path = temp_fasta.name
        
        # Create temporary output file for BLAST results
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_output:
            temp_output_path = temp_output.name
        
        # Run BLAST with hardcoded max_target_seqs=500 to get more hits for filtering
        blast_cmd = [
            'blastn',
            '-db', blast_db_path,
            '-query', temp_fasta_path,
            '-out', temp_output_path,
            '-evalue', str(evalue_threshold),
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle',
            '-max_target_seqs', '500',
            '-num_threads', str(threads)
        ]
        
        result = subprocess.run(blast_cmd, capture_output=True, text=True, check=True)
        
        # Parse BLAST results and filter to max_blast_hits
        region_blast_results = []
        with open(temp_output_path, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 12:
                        sseqid = fields[1]
                        # Parse the new FASTA header format: Order_Family_genus_Species_accession_feature_type_[gene_name]
                        # Extract hit_origin: everything before the feature_type_[gene_name] part (excluding feature_type)
                        hit_origin = '_'.join(sseqid.split('[')[0].split('_')[:-1])
                        hit_origin = '_'.join(hit_origin.split('_')[:-1])  # Remove feature_type
                        # Extract hit_gene: the gene name inside the brackets
                        hit_gene = sseqid.split('[')[-1].rstrip(']')
                                            
                        blast_result = {
                            'sample_name': sample_name,
                            'intergenic_region_id': fields[0],
                            'hit_origin': hit_origin,
                            'hit_gene': hit_gene,
                            'percent_identity': float(fields[2]),
                            'alignment_length': int(fields[3]),
                            'mismatches': int(fields[4]),
                            'gap_opens': int(fields[5]),
                            'query_start': int(fields[6]),
                            'query_end': int(fields[7]),
                            'subject_start': int(fields[8]),
                            'subject_end': int(fields[9]),
                            'evalue': float(fields[10]),
                            'bitscore': float(fields[11]),
                            'hit_title': fields[12] if len(fields) > 12 else '',
                            'intergenic_start': region['start'],
                            'intergenic_end': region['end'],
                            'intergenic_length': region['length'],
                            'upstream_gene': region['upstream_gene'],
                            'downstream_gene': region['downstream_gene']
                        }
                        region_blast_results.append(blast_result)
        
        # Filter to max_blast_hits (sorted by E-value, smallest first) and remove duplicates
        region_blast_results.sort(key=lambda x: x['evalue'])
        
        # Remove duplicates based on hit_origin and hit_gene combination
        seen_hits = set()
        unique_results = []
        for result in region_blast_results:
            hit_key = (result['hit_origin'], result['hit_gene'])
            if hit_key not in seen_hits:
                seen_hits.add(hit_key)
                unique_results.append(result)
                if len(unique_results) >= max_blast_hits:
                    break
        
        blast_results.extend(unique_results)    
        
        # Clean up temporary files
        os.unlink(temp_fasta_path)
        os.unlink(temp_output_path)
            
    return blast_results


def load_existing_blast_results(report_file, sample_name):
    """
    Load existing BLAST results from a TSV report file.
    
    Args:
        report_file (str): Path to the TSV report file
        sample_name (str): Name of the sample
        
    Returns:
        list: List of BLAST result dictionaries
    """
    blast_results = []
    
    try:
        with open(report_file, 'r') as f:
            # Skip header line
            next(f)
            
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 15:  # Ensure we have enough fields
                        blast_result = {
                            'sample_name': fields[0],
                            'intergenic_region_id': fields[1],
                            'hit_origin': fields[2],
                            'hit_gene': fields[3],
                            'percent_identity': float(fields[4]),
                            'alignment_length': int(fields[5]),
                            'mismatches': int(fields[6]),
                            'gap_opens': int(fields[7]),
                            'query_start': int(fields[8]),
                            'query_end': int(fields[9]),
                            'subject_start': int(fields[10]),
                            'subject_end': int(fields[11]),
                            'evalue': float(fields[12]),
                            'bitscore': float(fields[13]),
                            'hit_title': fields[14] if len(fields) > 14 else '',
                            'intergenic_start': int(fields[15]) if len(fields) > 15 else 0,
                            'intergenic_end': int(fields[16]) if len(fields) > 16 else 0,
                            'intergenic_length': int(fields[17]) if len(fields) > 17 else 0,
                            'upstream_gene': fields[18] if len(fields) > 18 else '',
                            'downstream_gene': fields[19] if len(fields) > 19 else ''
                        }
                        blast_results.append(blast_result)
    except Exception as e:
        # Use global logger if available, otherwise just raise the exception
        if 'logger' in globals():
            logger.error(f"{"[ERROR]:":10} Failed to load existing BLAST results from {report_file}: {str(e)}")
        raise
    
    return blast_results


def write_blast_report(blast_results, output_file, sample_name, logger):
    """
    Write BLAST results to TSV file.
    
    Args:
        blast_results (list): List of BLAST result dictionaries
        output_file (str): Output file path
        sample_name (str): Name of the sample or "ALL_GENOMES" for combined report
    """
    if not blast_results:
        logger.warning(f"{"[WARNING]:":10} No BLAST results to write for {sample_name}")
        raise Exception("No BLAST results to write")
    
    # Define column headers
    headers = [
        'sample_name', 'intergenic_region_id', 'hit_origin', 'hit_gene', 'percent_identity',
        'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end',
        'subject_start', 'subject_end', 'evalue', 'bitscore', 'hit_title',
        'intergenic_start', 'intergenic_end', 'intergenic_length', 'upstream_gene', 'downstream_gene'
    ]
    
    # Sort BLAST results by E-value (smallest first)
    sorted_blast_results = sorted(blast_results, key=lambda x: x.get('evalue', float('inf')))
    
    # Write TSV file
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(headers) + '\n')
        
        # Write data
        for result in sorted_blast_results:
            row = [str(result.get(header, '')) for header in headers]
            f.write('\t'.join(row) + '\n')
    
    logger.debug(f"{"[DEBUG]:":10} BLAST report written for {sample_name}: {output_file}")
        
 

def parse_metadata_tsv(metadata_file_path, input_file_dir, entry='main'):
    """
    Parse and validate the metadata TSV file for EMBL conversion.
    
    This function reads a TSV file containing sample metadata and validates:
    1. File format and required columns
    2. Data completeness for required fields (input_filename, linear_or_circular)
    3. Uses defaults for optional empty fields (project_id, locus_tag, genus_species)
    4. Validates that ALL INPUT files in the input directory have metadata entries if entry is not provided
    
    Args:
        metadata_file_path (str): Path to the metadata TSV file
        input_file_dir (str): Path to directory containing input files
        entry (str): Optional entry name to parse. 
        
    Returns:
        dict: Dictionary mapping input_filename to metadata dictionary containing:
            - project_id: Project identifier (uses default if empty)
            - locus_tag: Locus tag prefix (uses default if empty)
            - genus_species: Genus and species name (uses default if empty)
            - linear_or_circular: Whether the genome is linear or circular (required)
            
    Raises:
        SystemExit: If the file cannot be parsed or contains invalid data
    """

    logger.info(f"{"[INFO]:":10} Parsing metadata TSV file: {metadata_file_path}")

    required_columns = [
        'input_filename', 'project_id', 'locus_tag', 'genus_species', 'linear_or_circular'
    ]

    try:
        # Read the TSV file
        metadata_df = pd.read_csv(metadata_file_path, sep='\t', dtype=str)
        
        # Check for required columns
        if entry == 'main': # i.e. running main() for subcommand `annotate_and_check` - otherwise metadata is optional
            missing_columns = set(required_columns) - set(metadata_df.columns)
            if missing_columns:
                logger.error(f"{"[ERROR]:":10} Missing required columns in metadata TSV: {missing_columns}")
                logger.error(f"{"":10} Required columns: {required_columns}")
                utils.exit_program()
        
            # Check for empty dataframe
            if metadata_df.empty:
                logger.error(f"{"[ERROR]:":10} Metadata TSV file is empty")
                utils.exit_program()

            # Check for missing values in required fields (only input_filename and linear_or_circular must have values)
            required_fields = ['input_filename', 'linear_or_circular']
            for field in required_fields:
                missing_values = metadata_df[metadata_df[field].isna() | (metadata_df[field] == '')]
                if not missing_values.empty:
                    logger.error(f"{"[ERROR]:":10} Missing values in required field '{field}':")
                    for _, row in missing_values.iterrows():
                        logger.error(f"{"":10} {row['input_filename']}")
                    utils.exit_program()

            # Validate linear_or_circular values
            valid_topologies = {'linear', 'circular'}
            invalid_topologies = metadata_df[~metadata_df['linear_or_circular'].isin(valid_topologies)]
            if not invalid_topologies.empty:
                logger.error(f"{"[ERROR]:":10} Invalid values in 'linear_or_circular' field. Must be 'linear' or 'circular':")
                for _, row in invalid_topologies.iterrows():
                    logger.error(f"{"":10} input_filename: {row['input_filename']}, value: '{row['linear_or_circular']}'")
                utils.exit_program()
        
        # Check for duplicate input filenames
        duplicates = metadata_df[metadata_df.duplicated(['input_filename'], keep=False)]
        if not duplicates.empty:
            logger.error(f"{"[ERROR]:":10} Duplicate input_filename entries found:")
            for _, row in duplicates.iterrows():
                logger.error(f"{"":10} {row['input_filename']}")
            utils.exit_program()
        
        # Convert to dictionary, using defaults for empty optional fields
        metadata_dict = {}
        
        for _, row in metadata_df.iterrows():
            
            # Helper function to safely get and strip values, handling NaN
            def safe_get_strip(value, default):
                if pd.isna(value):
                    return default
                return str(value).strip() or default
            
            # Handle input_filename (required field)
            input_filename_raw = row['input_filename']
            if pd.isna(input_filename_raw):
                logger.warning(f"{"[WARNING]:":10} Skipping row with missing 'input_filename' field")
                continue
            input_filename = safe_get_strip(input_filename_raw, '')
            
            if not input_filename:
                logger.warning(f"{"[WARNING]:":10} Skipping row with empty 'input_filename' field")
                continue

            metadata_dict[input_filename] = {
                'project_id': safe_get_strip(row.get('project_id'), 'UNKNOWN_PROJECT'),
                'locus_tag': safe_get_strip(row.get('locus_tag'), 'DEFAULT_TAG'),
                'genus_species': safe_get_strip(row.get('genus_species'), 'Unknown species'),
                'linear_or_circular': safe_get_strip(row.get('linear_or_circular'), 'unknown')
            }
        
        logger.debug(f"{"[DEBUG]:":10} Successfully parsed metadata TSV file: {len(metadata_dict)} samples")
        logger.debug(f"{"[DEBUG]:":10} Sample entries: {list(metadata_dict.keys())[:5]}")
        
        # Check that all INPUT files in the input directory have metadata entries
        if entry == 'main':
            input_files = glob.glob(os.path.join(input_file_dir, "*.fasta")) + \
                     glob.glob(os.path.join(input_file_dir, "*.fa")) + \
                     glob.glob(os.path.join(input_file_dir, "*.fas"))
            
        elif entry == 'check':
            input_files = glob.glob(os.path.join(input_file_dir, "*.gb")) + \
                     glob.glob(os.path.join(input_file_dir, "*.gbk")) + \
                     glob.glob(os.path.join(input_file_dir, "*.gb.gz")) + \
                     glob.glob(os.path.join(input_file_dir, "*.gbk.gz"))
        else:
            input_files = glob.glob(os.path.join(input_file_dir, "*"))  # not specific to file extension
        
        if not input_files:
            logger.warning(f"{"[WARNING]:":10} No INPUT files found in {input_file_dir}")
        else:
            # Get list of INPUT filenames that should be in metadata
            expected_input_files = []
            for input_file in input_files:
                input_filename = os.path.basename(input_file)
                expected_input_files.append(input_filename)
            
            # Check for missing metadata entries
            missing_from_metadata = []
            for input_filename in expected_input_files:
                if input_filename not in metadata_dict:
                    missing_from_metadata.append(input_filename)
            
            logger.info(f"{"[INFO]:":10} Metadata coverage: "
                        f"{len(expected_input_files) - len(missing_from_metadata)}/{len(expected_input_files)} samples")
            
            if entry == 'main': # i.e. running main() for subcommand `annotate_and_check` - otherwise metadata is optional
                if missing_from_metadata:
                    logger.error(f"{"[ERROR]:":10} ALL samples must be present in metadata file!")
                    logger.error(f"{"":10} {len(missing_from_metadata)} samples are missing from metadata:")
                    for input_filename in missing_from_metadata:
                        logger.error(f"{"":15} {input_filename}")
                    logger.error(f"{"":10} Please add all missing samples to your metadata TSV file")
                    utils.exit_program()
        
        utils.log_separator(logger)
        return metadata_dict
        
    except pd.errors.EmptyDataError:
        logger.error(f"{"[ERROR]:":10} Metadata TSV file is empty or contains no data")
        utils.exit_program()
    except pd.errors.ParserError as e:
        logger.error(f"{"[ERROR]:":10} Error parsing metadata TSV file: {e}")
        utils.exit_program()
    except Exception as e:
        logger.error(f"{"[ERROR]:":10} Unexpected error reading metadata TSV file: {e}")
        logger.error(traceback.format_exc())
        utils.exit_program()


def load_annotated_genbank_files(genbank_dir, output_directory, logger):
    """
    Load already annotated GenBank files into the annotated_genomes_dict structure.
    
    Args:
        genbank_dir (str): Directory containing annotated GenBank files
        output_directory (str): Output directory for split GenBank files
        logger: Logger instance for logging messages
        
    Returns:
        dict: Dictionary with sample names as keys and dicts containing:
            - 'gbk': Path to GenBank file (single) or list of paths (multi-sequence)
            - 'input_filename': Original input filename (placeholder)
            - 'is_multi_sequence': Boolean indicating if this is a multi-sequence sample
            - 'sequence_count': Number of sequences in this sample
              matching the structure expected by check_genes()
    """ 
    
    annotated_genomes_dict = {}
    
    # Check if directory exists
    if not os.path.exists(genbank_dir):
        logger.error(f"{"[ERROR]:":10} GenBank directory does not exist: {genbank_dir}")
        utils.exit_program()
    
    # Create output directory for individual files written from multi-record GenBank files
    output_dir = os.path.join(output_directory, "01_gbk_split")
    os.makedirs(output_dir, exist_ok=True)
    logger.debug(f"{"[DEBUG]:":10} Created output directory: {output_dir}")
    
    # Get all .gb and .gbk files in the directory
    genbank_files = []
    for ext in ['*.gb', '*.gbk', '*.gb.gz', '*.gbk.gz']:
        genbank_files.extend(glob.glob(os.path.join(genbank_dir, ext)))
    
    if not genbank_files:
        logger.error(f"{"[ERROR]:":10} No GenBank files (.gb, .gbk, .gb.gz, or .gbk.gz) found in directory: {genbank_dir}")
        utils.exit_program()
    
    logger.info(f"{"[INFO]:":10} Found {len(genbank_files)} GenBank files")
    
    # Process each GenBank file and split multi-record files
    sample_files = {}
    for gbk_file in genbank_files:
        try:
            # Extract sample name from filename (remove .gz and GenBank extensions)
            base_filename = os.path.basename(gbk_file)
            if base_filename.endswith('.gz'):
                # Remove .gz first, then remove GenBank extension
                sample_name = os.path.splitext(os.path.splitext(base_filename)[0])[0]
            else:
                # Just remove GenBank extension
                sample_name = os.path.splitext(base_filename)[0]
            
            # Initialize sample in dictionary if not present
            if sample_name not in sample_files:
                sample_files[sample_name] = {'gbk': [],
                                             'input_filename': base_filename}
            
            # Parse the GenBank file (may contain multiple records)
            # Parse GenBank file using utility function with warning filtering
            seq_records = utils.parse_genbank_file(gbk_file, logger)

            if not seq_records:
                logger.warning(f"{"[WARNING]:":10} No sequences found in {gbk_file}")
                continue
            
            # If single record, use original file
            if len(seq_records) == 1:
                sample_files[sample_name]['gbk'].append(gbk_file)
            
            else:
                # Multiple records - split into individual files
                logger.info(f"{"[INFO]:":10} Splitting {len(seq_records)} records from {base_filename}")

                for i, record in enumerate(seq_records, 1):
                    # Create individual GenBank file
                    # Create a safe sequence name (replace spaces and special chars)
                    safe_seq_name = str(record.id).replace(' ', '_').replace(':', '_').replace('|', '_')
                    individual_gbk_name = f"{sample_name}_seq{i:03d}_{safe_seq_name}.gb"
                    individual_gbk_path = os.path.join(output_dir, individual_gbk_name)
                    
                    # Write individual GenBank record
                    SeqIO.write(record, individual_gbk_path, 'genbank')
                    
                    # Add to sample files
                    sample_files[sample_name]['gbk'].append(individual_gbk_path)
                    logger.debug(f"{"[DEBUG]:":10} Created {individual_gbk_name}")
                    
        except Exception as e:
            logger.error(f"{"[ERROR]:":10} Error processing GenBank file {gbk_file}: {e}")
            utils.log_manager.handle_error(e, traceback.format_exc(), "load_annotated_genbank_files()")
    
    # Convert to the expected structure
    for sample_name, files in sample_files.items():
        gbk_files = files['gbk']
        input_filename = files['input_filename']
        
        # Determine if this is a multi-sequence sample
        is_multi_sequence = len(gbk_files) > 1
        sequence_count = len(gbk_files)
        
        # Create the sample data structure matching annotate_genomes() output
        sample_data = {
            'input_filename': input_filename,  
            'is_multi_sequence': is_multi_sequence,
            'sequence_count': sequence_count
        }
        
        # Add GenBank file paths (convert to single path if only one file, keep as list if multiple)
        if len(gbk_files) == 1:
            sample_data['gbk'] = gbk_files[0]
        else:
            sample_data['gbk'] = gbk_files
        
        annotated_genomes_dict[sample_name] = sample_data
        logger.debug(f"{"[DEBUG]:":10} Loaded {sequence_count} sequences for {sample_name}")

    logger.info(f"{"[INFO]:":10} Successfully loaded {len(annotated_genomes_dict)} samples with annotated GenBank files")

    utils.log_separator(logger)
    time.sleep(0.1)
    
    return annotated_genomes_dict


def check_pipeline(args):
    """Continue pipeline from annotated GenBank files.
    Args:
        args (argparse.Namespace): Parsed command line arguments containing input file paths
            and other configuration options.

    Returns:
        None: No return value specified.

    Raises:
        SystemExit: If input files or directories do not exist.
    """

    # Track wall-clock runtime for completion message
    start_time = time.time()

    try:
        global logger, log_queue, log_listener

        # Set up log and report directories
        utils.setup_log_and_report_directories(args)

        # Set up global logger
        logger, log_queue, log_listener = utils.log_manager.setup(
            __name__, 'check', log_directory=args.log_directory
        )

        # Print arguments to screen and log:
        utils.print_arguments(args, logger, __version__)

        # Check for external dependencies:
        utils.check_dependencies(args, logger, entry='check')

        # Resolve base data directory after dependency check
        data_dir_base = utils.resolve_data_dir_base()

        # Load gene median lengths from resolved data directory
        gene_median_lengths = load_gene_median_lengths(data_dir_base=data_dir_base)

        # Load gene synonyms from resolved data directory
        gene_synonyms = load_gene_synonyms(data_dir_base=data_dir_base)

        # Check no_alignment and refs_order
        ref_gene_seqrecords = get_references(args, gene_median_lengths, gene_synonyms, data_dir_base=data_dir_base)

        # Parse required metadata TSV file
        if args.metadata_tsv:
            metadata_dict = parse_metadata_tsv(args.metadata_tsv, args.annotated_genbank_dir, entry='check')
        else:
            metadata_dict = None
            logger.info(f"{"[INFO]:":10} No metadata TSV file provided, default values with be used in EMBL output")
            utils.log_separator(logger)
            time.sleep(0.1)

        # Load annotated GenBank files
        annotated_genomes_dict = load_annotated_genbank_files(args.annotated_genbank_dir, args.output_directory, logger)

        # Check genes and write reports
        all_sample_results = check_genes(gene_median_lengths, annotated_genomes_dict, args.min_length_percentage,
                                         args.max_length_percentage, args.report_directory, log_queue, args.pool,
                                         gene_synonyms)

        # Convert assembly gbk files to embl format
        convert_gbk_to_embl(annotated_genomes_dict, args.output_directory, metadata_dict=metadata_dict)

        # Generate alignments if not disabled
        if not args.no_alignment:
            align_genes(all_sample_results, ref_gene_seqrecords, args.output_directory, args.pool, args.threads,
                        args.refs_order)

        # Query intergenic regions
        if not args.skip_intergenic_analysis:
            query_intergenic_regions(annotated_genomes_dict, args.output_directory, args.min_intergenic_length,
                                     args.blast_evalue, args.debug_intergenic, args.max_blast_hits, args.pool,
                                     args.threads, log_queue, args.custom_blast_db)
        else:
            logger.info(f"{"[INFO]:":10} Skipping intergenic region analysis as requested")
 
    except Exception as e:
        utils.log_manager.handle_error(e, traceback.format_exc(), "check_pipeline()")

    finally:
        # Log total completion time before cleaning up the logger
        utils.log_separator(logger)
        utils.log_completion_time(start_time, logger if ('logger' in globals() and logger) else None,
                                  label="PAV subcommand `check` completed")

        utils.log_manager.cleanup()


def main(args):
    """Main annotation and validation pipeline for plastid genomes.
    
    This function orchestrates the complete PAV pipeline including:
    1. Genome annotation using Chloë
    2. Gene validation against reference data
    3. EMBL format conversion
    4. Multiple sequence alignment generation
    5. Intergenic region analysis
    
    Args:
        args (argparse.Namespace): Parsed command line arguments containing input file paths
            and other configuration options.

    Returns:
        None: No return value specified.

    Raises:
        SystemExit: If input files or directories do not exist.
    """

    # Track wall-clock runtime for completion message
    start_time = time.time()

    try:
        global logger, log_queue, log_listener

        # Set up log and report directories
        utils.setup_log_and_report_directories(args)

        # Set up global logger
        logger, log_queue, log_listener = utils.log_manager.setup(
            __name__, 'annotate_and_check', log_directory=args.log_directory
        ) 

        # Print arguments to screen and log:
        utils.print_arguments(args, logger, __version__)

        # Check for external dependencies:
        utils.check_dependencies(args, logger)

        # Resolve base data directory after dependency check
        data_dir_base = utils.resolve_data_dir_base()

        # Load gene median lengths from resolved data directory
        gene_median_lengths = load_gene_median_lengths(data_dir_base=data_dir_base)

        # Load gene synonyms from resolved data directory
        gene_synonyms = load_gene_synonyms(data_dir_base=data_dir_base)

        # Check no_alignment and refs_order
        ref_gene_seqrecords = get_references(args, gene_median_lengths, gene_synonyms, data_dir_base=data_dir_base)

        # Parse required metadata TSV file
        metadata_dict = parse_metadata_tsv(args.metadata_tsv, args.genome_fasta_dir)

        # Validate linearisation genes against gene_synonyms.txt
        validated_linearise_genes = validate_linearisation_genes(args.linearise_gene, gene_synonyms, logger)

        # Annotate the genomes using chloë (project dir required)
        annotated_genomes_dict = annotate_genomes(
            args.genome_fasta_dir,
            args.output_directory,
            args.chloe_project_dir,
            validated_linearise_genes,
            metadata_dict,
            args.pool,
            log_queue
        )
        
        # Check genes and write reports
        all_sample_results = check_genes(gene_median_lengths, annotated_genomes_dict, args.min_length_percentage,
                                         args.max_length_percentage, args.report_directory, log_queue, args.pool,
                                         gene_synonyms)

        # Convert assembly gbk files to embl format
        convert_gbk_to_embl(annotated_genomes_dict, args.output_directory, metadata_dict=metadata_dict)

        # Generate alignments if not disabled
        if not args.no_alignment:
            align_genes(all_sample_results, ref_gene_seqrecords, args.output_directory, args.pool, args.threads,
                        args.refs_order)

        # Query intergenic regions
        if not args.skip_intergenic_analysis:
            query_intergenic_regions(annotated_genomes_dict, args.output_directory, args.min_intergenic_length,
                                     args.blast_evalue, args.debug_intergenic, args.max_blast_hits, args.pool,
                                     args.threads, log_queue, args.custom_blast_db)
        else:
            logger.info(f"{"[INFO]:":10} Skipping intergenic region analysis as requested")
 
    except Exception as e:
        utils.log_manager.handle_error(e, traceback.format_exc(), "main()")

    finally:
        # Log total completion time before cleaning up the logger
        utils.log_separator(logger)
        utils.log_completion_time(start_time, logger if ('logger' in globals() and logger) else None,
                                  label="PAV subcommand `annotate_and_check` completed")

        utils.log_manager.cleanup()
