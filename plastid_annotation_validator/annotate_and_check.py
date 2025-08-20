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
import gzip
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
import traceback
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import subprocess
import tempfile
import glob
import re
import shutil
import uuid
from Bio import AlignIO
import io
import pandas as pd
import importlib.resources
import time
from pathlib import Path

# Import NewTargets modules:
from plastid_annotation_validator.version import __version__
from plastid_annotation_validator import utils

# Initialise logger objects to None
logger = None
log_queue = None
log_listener = None

def load_gene_synonyms():
    """
    Load gene synonyms from the package text file.
    
    Args:
        None

    Returns:
        dict: Dictionary mapping gene names to their standardized synonyms
        
    """

    with importlib.resources.open_text('plastid_annotation_validator.data', 'gene_synonyms.txt') as synonyms_file:
        gene_synonyms = {}
        
        for line in synonyms_file:
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Parse the mapping: "GenBank_name -> Standard_name"
            if ' = ' in line:
                genbank_name, standard_name = line.split(' = ', 1)
                gene_synonyms[genbank_name.strip()] = standard_name.strip()

    logger.debug(f"{"[DEBUG]:":10} Successfully loaded {len(gene_synonyms)} gene synonyms from package resources")
    return gene_synonyms
        

def load_gene_median_lengths():
    """ 
    Load gene median lengths from the package CSV file.
    
    Args:
        logger: Logger instance for logging messages

    Returns:
        dict: Dictionary mapping gene names to their median lengths
        
    """

    with importlib.resources.open_text('plastid_annotation_validator.data', 'plDNA_genes_median_lengths.csv') as csv_file:
        gene_lengths_df = pd.read_csv(csv_file)

    logger.debug(f"{"[DEBUG]:":10} Successfully loaded reference gene median lengths from package resources")
    logger.debug(f"{"[DEBUG]:":10} Loaded {len(gene_lengths_df)} gene entries")
    
    # Convert to dictionary for easier access
    gene_median_lengths = dict(zip(gene_lengths_df['Key'], gene_lengths_df['Value']))
    
    logger.debug(f"{"[DEBUG]:":10} Loaded {len(gene_median_lengths)} gene median lengths")
    
    return gene_median_lengths


def parse_gbk_genes(gbk_file_path, fasta_file_path, logger, gene_synonyms=None):
    """
    Parse gene annotations from a GenBank file and extract gene lengths, CDS, rRNA, and tRNA information.
    
    This function processes a GenBank file to extract comprehensive gene annotation information:
    1. Parses CDS, rRNA, and tRNA features from the GenBank file
    2. Extracts gene names and applies synonym mapping for standardization
    3. Calculates gene lengths and handles multiple gene copies
    4. Extracts sequence information for translation validation
    5. Organizes information into structured dictionaries for downstream analysis
    
    Args:
        gbk_file_path (str): Path to the GenBank file (.gb or .gbk)
        fasta_file_path (str): Path to the corresponding FASTA file for sequence length validation
        logger: Logger instance for logging messages and debugging
        gene_synonyms (dict, optional): Dictionary mapping gene names to standardized synonyms.
            Used to normalize gene naming across different annotation sources.
        
    Returns:
        tuple: (gene_lengths, gene_cds_info, gene_rRNA_info, gene_tRNA_info) where:
            - gene_lengths (dict): Dictionary mapping gene names to length information containing:
                - 'length': Representative gene length (mean if multiple copies)
                - 'copies': Number of copies of this gene in the genome
                - 'mean_length': Mean length across all copies
                - 'copy_lengths': List of individual copy lengths
            - gene_cds_info (dict): Dictionary mapping gene names to CDS information for translation checking
            - gene_rRNA_info (dict): Dictionary mapping gene names to rRNA information
            - gene_tRNA_info (dict): Dictionary mapping gene names to tRNA information
    
    Raises:
        ValueError: If no records are found in the GenBank file or if a feature lacks a gene name
        Exception: Various exceptions that may occur during file parsing or sequence extraction
        
    Note:
        The function handles multiple copies of the same gene by calculating mean lengths
        and storing individual copy information. Gene names are standardized using the
        provided synonym mapping if available.
    """
    gene_lengths = {}
    gene_cds_info = {}
    gene_rRNA_info = {}
    gene_tRNA_info = {}
    
    try:
        # Create adjusted GenBank file object with corrected LOCUS line
        gbk_io, _ = utils.build_adjusted_genbank_io(
            gbk_file_path=gbk_file_path,
            fasta_file_path=fasta_file_path,
            logger=logger,
        )
        
        # Parse GenBank file using Biopython
        records = list(SeqIO.parse(gbk_io, "genbank"))
        if not records:
            raise ValueError(f"No records found in GenBank file: {gbk_file_path}")
        
        # Process each record (usually only one for plastid genomes)
        for record in records:
            # Initialize data structures for this record
            gene_features = defaultdict(list)
            gene_cds_locations = defaultdict(list)
            gene_rRNA_locations = defaultdict(list)
            gene_tRNA_locations = defaultdict(list)
            
            # Process each feature in the record
            for feature in record.features:
                if feature.type not in ["CDS", "rRNA", "tRNA"]:
                    continue
                    
                # Extract gene name from qualifiers
                gene_name = _extract_gene_name(feature)
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
                    logger.debug(f"{"[DEBUG]:":10} Gene {mapped_gene_name} length mismatch: Expected {gene_length} != Extracted {len(extracted_seq)}")
                    length_mismatch_msg = f"Gene {mapped_gene_name} length mismatch: Expected {gene_length} != Extracted {len(extracted_seq)}"
                
                # Store gene length information
                gene_features[mapped_gene_name].append(gene_length)
                
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
            
            # Calculate final gene lengths and organize information
            for gene_name, lengths in gene_features.items():
                # Calculate representative length (mean if multiple copies)
                if len(lengths) == 1:
                    representative_length = lengths[0]
                else:
                    representative_length = round(sum(lengths) / len(lengths))
                    logger.debug(f"{"[DEBUG]:":10} Gene {gene_name} has {len(lengths)} copies with lengths: {lengths}, using mean: {representative_length:.1f}")
                
                # Store gene length information
                gene_lengths[gene_name] = {
                    'length': representative_length,
                    'copies': len(lengths),
                    'mean_length': representative_length,
                    'copy_lengths': lengths
                }
                
                # Store CDS information if available
                if gene_name in gene_cds_locations:
                    gene_cds_info[gene_name] = gene_cds_locations[gene_name]
                
                # Store rRNA information if available
                if gene_name in gene_rRNA_locations:
                    gene_rRNA_info[gene_name] = gene_rRNA_locations[gene_name]
                
                # Store tRNA information if available
                if gene_name in gene_tRNA_locations:
                    gene_tRNA_info[gene_name] = gene_tRNA_locations[gene_name]
        
        logger.debug(f"{"[DEBUG]:":10} Parsed {len(gene_lengths)} genes (CDS, rRNA, and tRNA) from GenBank file")
        return gene_lengths, gene_cds_info, gene_rRNA_info, gene_tRNA_info
        
    except Exception as e:
        logger.error(f"{"[ERROR]:":10} Error parsing GenBank file: {e}")
        logger.error(traceback.format_exc())
        raise


def _extract_gene_name(feature):
    """
    Extract gene name from a Biopython feature's qualifiers.
    
    Args:
        feature: Biopython SeqFeature object
        
    Returns:
        str: Gene name if found, None otherwise
    """
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    elif 'product' in feature.qualifiers:
        # Use first word of product name as fallback
        return feature.qualifiers['product'][0].split()[0]
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
                warnings.append(f"Start codon is {start_codon}, expected ATG for copy {i+1}")
            
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


def check_single_sample_genes(sample_data, gene_median_lengths, min_threshold, max_threshold, gene_synonyms=None, log_queue=None):
    """
    Validate gene annotations for a single plastid genome sample by checking gene lengths and translation quality.
    
    This function performs comprehensive validation of gene annotations for one sample by:
    1. Parsing the GenBank file to extract gene annotations and sequences
    2. Comparing gene lengths against reference median lengths from a curated database
    3. Checking translation quality for CDS genes to identify potential annotation issues
    4. Applying gene synonym mapping for standardized gene naming
    5. Generating detailed validation results for each gene
    
    Args:
        sample_data (dict): Dictionary containing sample information with keys:
            - 'sample_name': Name of the sample being processed
            - 'gbk_file': Path to the annotated GenBank file
            - 'fasta_file': Path to the corresponding FASTA file
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
                - 'gene_results': Dictionary mapping gene names to validation results containing:
                    - 'gene_name': Gene identifier
                    - 'actual_length': Observed gene length in base pairs
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
                - 'gene_cds_info': CDS location information for translation checking
                - 'gene_rRNA_info': rRNA location information
                - 'gene_tRNA_info': tRNA location information
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
        gbk_file = sample_data['gbk_file']
        fasta_file = sample_data['fasta_file']

        length_warnings = set()
        sample_dict = {}
    
        # Parse gene lengths and CDS information from GenBank file (adjusted in-memory)
        gene_lengths, gene_cds_info, gene_rRNA_info, gene_tRNA_info = parse_gbk_genes(gbk_file, fasta_file, worker_logger, gene_synonyms)

        # Check for genes in CDS info that are not in gene lengths
        cds_only_genes = set(gene_cds_info.keys()) - set(gene_lengths.keys())
        if cds_only_genes:
            return False, (ValueError(f"Genes in CDS info but not in gene lengths: {list(cds_only_genes)[:10]}"), traceback.format_exc())
        
        # Check for genes present in gene_lengths but not in gene_median_lengths, and vice versa
        genes_in_lengths_not_in_median = list(set(gene_lengths.keys()) - set(gene_median_lengths.keys()))
        genes_in_median_not_in_lengths = list(set(gene_median_lengths.keys()) - set(gene_lengths.keys()))
        
        if genes_in_lengths_not_in_median:
            worker_logger.debug(f"{"[DEBUG]:":10} Genes in sample {sample_name} but not in median lengths: {genes_in_lengths_not_in_median}...")
            # tqdm.write(f"{"[INFO]:":10} Genes in sample {sample_name} but not in median lengths: {genes_in_lengths_not_in_median}...")
        if genes_in_median_not_in_lengths:
            worker_logger.debug(f"{"[DEBUG]:":10} Genes in median lengths but not in sample {sample_name}: {genes_in_median_not_in_lengths}...")
            # tqdm.write(f"{"[INFO]:":10} Genes in median lengths but not in sample {sample_name}: {genes_in_median_not_in_lengths}...")
       
        # Check each gene against median lengths and translation
        for gene_name, gene_info in gene_lengths.items():

            actual_length = gene_info['length']
            copies = gene_info['copies']
            copy_lengths = gene_info['copy_lengths']

            # Add gene to sample_dict with basic info and placeholders
            sample_dict[gene_name] = {
                'gene_name': gene_name,
                'actual_length': actual_length,
                'copies': copies,
                'copy_lengths': copy_lengths,
                'threshold_min': min_threshold,
                'threshold_max': max_threshold,
                'median_length': None,
                'issue': None,
                'details': None,
                'translation_status': 'N/A',
                'translation_details': 'N/A',
            }
            
            # Check length against median
            if gene_name in gene_median_lengths:
                median_length = gene_median_lengths[gene_name]
                min_expected = median_length * min_threshold
                max_expected = median_length * max_threshold

                # Add median length info to sample_dict
                sample_dict[gene_name]['median_length'] = median_length
   
                # Check if gene is too short
                if actual_length < min_expected:
                    warning_msg = f"Too short: {actual_length} bp < {min_expected:.0f} bp ({min_threshold*100}% of ref median)"
                    sample_dict[gene_name]['issue'] = 'too_short'
                    sample_dict[gene_name]['details'] = warning_msg
                    length_warnings.add(gene_name)

                # Check if gene is too long
                elif actual_length > max_expected:
                    warning_msg = f"Too long: {actual_length} bp > {max_expected:.0f} bp ({max_threshold*100}% of median)"
                    sample_dict[gene_name]['issue'] = 'too_long'
                    sample_dict[gene_name]['details'] = warning_msg
                    length_warnings.add(gene_name)

                else:
                    sample_dict[gene_name]['issue'] = 'OK'
                    sample_dict[gene_name]['details'] = 'OK'

            else:
                worker_logger.warning(f"{"[WARNING]:":10} Gene {gene_name} not in gene_median_lengths")
                # tqdm.write(f"{"[WARNING]:":10} Gene {gene_name} not in gene_median_lengths")
                sample_dict[gene_name]['issue'] = 'not_in_median_lengths'
                sample_dict[gene_name]['details'] = 'Not in median lengths'
                length_warnings.add(gene_name)
            
            # Check translation for CDS genes and create alignments with reference CDSs (optional)
            if gene_name in gene_cds_info:
                translation_issues = check_gene_translation(gene_name, gene_cds_info[gene_name], logger)

                if translation_issues:
                    sample_dict[gene_name]['translation_status'] = 'FAIL'
                    sample_dict[gene_name]['translation_details'] = translation_issues
                    length_warnings.add(gene_name)
                else:
                    sample_dict[gene_name]['translation_status'] = 'OK'
                    sample_dict[gene_name]['translation_details'] = 'OK'

        # Collate results for return
        result_dict = {
            'sample_name': sample_name,
            'gene_results': sample_dict,
            'total_genes': len(gene_lengths),
            'genes_with_warnings': len(length_warnings),
            'gene_cds_info': gene_cds_info,
            'gene_rRNA_info': gene_rRNA_info,
            'gene_tRNA_info': gene_tRNA_info,
            'missing_genes': genes_in_median_not_in_lengths
        }

        return True, result_dict
        
    except Exception as e:
        return False, (e, traceback.format_exc())


def check_genes(gene_median_lengths, annotated_genomes_dict, min_threshold, max_threshold, 
                report_directory, log_queue, pool_size=1,  gene_synonyms=None):
    """
    Validate gene annotations across multiple plastid genomes by comparing gene lengths against reference median values.
    
    This function processes annotated plastid genomes and validates gene annotations by:
    1. Parsing GenBank files to extract gene annotations and sequences
    2. Comparing gene lengths against reference median lengths from a curated database
    3. Applying gene synonym mapping for standardized gene naming
    4. Generating comprehensive reports of validation results
    5. Using multiprocessing for efficient parallel processing
    
    Args:
        gene_median_lengths (dict): Dictionary mapping gene names to their reference median lengths
        annotated_genomes_dict (dict): Dictionary mapping sample names to dicts containing:
            - 'gbk': Path to annotated GenBank file
            - 'fasta': Path to corresponding FASTA file
        min_threshold (float): Minimum acceptable length as percentage of median (e.g., 80.0 for 80%)
        max_threshold (float): Maximum acceptable length as percentage of median (e.g., 120.0 for 120%)
        report_directory (str): Directory path where validation reports will be written
        log_queue (queue.Queue): Multiprocessing-safe queue for logging messages
        pool_size (int, optional): Number of worker processes for parallel processing. Defaults to 1.
        gene_synonyms (dict, optional): Dictionary mapping gene names to standardized synonyms.
            Used to normalize gene naming across different annotation sources.

    Returns:
        dict: Dictionary mapping sample names to validation results containing:
            - 'gene_results': Dictionary of gene validation results
            - 'genes_with_warnings': Count of genes that failed validation
            - 'total_genes': Total number of genes processed
            - 'error': Error message if processing failed (optional)
    
    Raises:
        SystemExit: If any sample fails to process and error handling is triggered
        
    Note:
        The function uses multiprocessing to process samples in parallel, with each sample
        being validated independently. Results are collected and combined into comprehensive
        reports including both individual sample reports and a combined summary.
    """

    all_sample_results = {}
    total_warnings = 0
    
    # Prepare data for multiprocessing
    sample_data_list = []
    for sample_name, data in annotated_genomes_dict.items():
        sample_data = {
            'fasta_file': data['fasta'],
            'gbk_file': data['gbk'],
            'sample_name': sample_name
        }
        sample_data_list.append(sample_data)
    
    # Process samples with multiprocessing and progress bar
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
        for future in tqdm(as_completed(future_to_sample), total=len(future_to_sample), desc=f"{"[INFO]:":10} {"Validating genes":<20}", file=sys.stdout):
            sample_name = future_to_sample[future]

            success, result = future.result()
            
            if success:

                all_sample_results[sample_name] = result
                total_warnings += result['genes_with_warnings']

            else:
                utils.log_manager.handle_error(result[0], "Gene length checking", sample_name)

    # Generate and write report
    write_gene_length_report(all_sample_results, logger, min_threshold, max_threshold, report_directory, gene_median_lengths)
    
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
    
    The reports include information about:
    - Gene length validation against reference median values
    - Translation quality checks for CDS genes
    - Multiple gene copy handling and analysis
    - Warning categorization and detailed descriptions
    
    Args:
        all_results (dict): Dictionary mapping sample names to validation results containing:
            - 'gene_results': Dictionary of gene validation results for each gene
            - 'error': Error message if sample processing failed (optional)
        logger: Logger instance for logging messages and progress updates
        min_threshold (float): Minimum acceptable length as percentage of median used for validation
        max_threshold (float): Maximum acceptable length as percentage of median used for validation
        report_directory (str): Directory path where all reports will be written
        gene_median_lengths (dict): Dictionary mapping gene names to their reference median lengths
    
    Returns:
        None: Reports are written directly to files in the specified directory
        
    Generated Files:
        - Individual sample reports: "{sample_name}_gene_validation_report.tsv"
        - Combined report: "all_samples_gene_validation_report.tsv"
        
    Report Format (TSV with headers):
        - Sample: Sample identifier
        - Gene_ID: Gene identifier
        - Gene_Name: Gene name (same as Gene_ID)
        - Actual_Length: Observed gene length in base pairs
        - Median_Length_Refs: Reference median length (NA if not available)
        - Min_Threshold: Minimum threshold percentage used
        - Max_Threshold: Maximum threshold percentage used
        - Copies: Number of copies of this gene in the genome
        - Copy_Lengths: Comma-separated list of individual copy lengths
        - Warning: Validation status and warning messages
        - Translation_Details: Details about translation quality issues
        
    Warning Categories:
        - 'too_short': Gene length below minimum threshold
        - 'too_long': Gene length above maximum threshold
        - 'not_in_median_lengths': Gene not found in reference database
        - 'Translation issue': Problems with CDS translation
        - 'OK': No validation issues found
        
    Note:
        The function handles both successful samples and samples with processing errors.
        Error samples are included in the combined report with appropriate error indicators.
        All reports use tab-separated values (TSV) format for easy parsing and analysis.
    """

    logger.info(f"{"[INFO]:":10} Writing gene validation reports to: {report_directory}\n")
    
    # Prepare combined TSV data
    combined_tsv_lines = [
        "Sample\tGene_ID\tGene_Name\tActual_Length\tMedian_Length_Refs\tMin_Threshold\tMax_Threshold\tCopies\tCopy_Lengths\tWarning\tTranslation_Details"]

    total_warnings_in_reports = 0

    # Process each sample
    for sample_name, result in all_results.items():
        if 'error' in result:
            # Add error row to combined report
            combined_tsv_lines.append(f"{sample_name}\tERROR\tERROR\t0\t0\t{min_threshold}\t{max_threshold}\t0\t\t{result['error']}\tERROR")
            continue
        
        # Get gene results from the result (this is the correct structure)
        gene_results = result.get('gene_results', {})
        
        logger.debug(f"{"[DEBUG]:":10} Processing {sample_name}: found {len(gene_results)} genes")
        
        # Prepare individual sample report
        sample_tsv_lines = [
            "Gene_ID\tGene_Name\tActual_Length\tMedian_Length_Refs\tMin_Threshold\tMax_Threshold\tCopies\tCopy_Lengths\tWarning\tTranslation_Details"]

        # Process all genes for this sample
        logger.debug(f"{"[DEBUG]:":10} Processing {len(gene_results)} genes for {sample_name}")
        warnings_found = 0
        missing_genes_count = 0
        
        for gene_name, gene_info in sorted(gene_results.items()):
            actual_length = gene_info['actual_length']
            copies = gene_info['copies']
            copy_lengths = gene_info['copy_lengths']
            copy_lengths_str = ','.join(map(str, copy_lengths))
            
            # Get median length and determine warning status
            median_length = gene_info.get('median_length')
            issue = gene_info.get('issue', 'unknown')
            details = gene_info.get('details', '')
            translation_status = gene_info.get('translation_status', 'OK')
            translation_details = gene_info.get('translation_details', '')
            
            # Determine warning message based on issue and translation status

            warning_msg = []
            if issue == 'too_short':
                warning_msg.append(details)
                warnings_found += 1
                total_warnings_in_reports += 1
            if issue == 'too_long':
                warning_msg.append(details)
                warnings_found += 1
                total_warnings_in_reports += 1
            if issue == 'not_in_median_lengths':
                warning_msg.append(details)
                warnings_found += 1
                total_warnings_in_reports += 1
            if translation_status == 'FAIL':
                warning_msg.append("Translation issue")
                warnings_found += 1
                total_warnings_in_reports += 1

            if warning_msg:
                warning_msg = ", ".join(warning_msg)
            else:
                warning_msg = "OK"
            
            # Add to both combined and individual reports
            combined_tsv_lines.append(f"{sample_name}\t{gene_name}\t{gene_name}\t{actual_length}\t{median_length or 'NA'}\t{min_threshold}\t{max_threshold}\t{copies}\t{copy_lengths_str}\t{warning_msg}\t{translation_details}")
            sample_tsv_lines.append(f"{gene_name}\t{gene_name}\t{actual_length}\t{median_length or 'NA'}\t{min_threshold}\t{max_threshold}\t{copies}\t{copy_lengths_str}\t{warning_msg}\t{translation_details}")
        
        # Add missing genes to reports
        missing_genes = result.get('missing_genes', [])
        for missing_gene in missing_genes:
            median_length = gene_median_lengths.get(missing_gene, 'NA')
            # Add missing gene to both combined and individual reports
            combined_tsv_lines.append(f"{sample_name}\t{missing_gene}\t{missing_gene}\tNA\t{median_length}\t{min_threshold}\t{max_threshold}\tNA\tNA\tMissing\tNA")
            sample_tsv_lines.append(f"{missing_gene}\t{missing_gene}\tNA\t{median_length}\t{min_threshold}\t{max_threshold}\tNA\tNA\tMissing\tNA")
            missing_genes_count += 1
            total_warnings_in_reports += 1
        
        # Write individual sample report
        sample_report_file = os.path.join(report_directory, f"{sample_name}_gene_validation_report.tsv")
        with open(sample_report_file, 'w') as f:
            f.write("\n".join(sample_tsv_lines))
        
        logger.info(f"{" ":15} Sample {sample_name}: ({len(sample_tsv_lines)-1} genes, {warnings_found} warnings, {missing_genes_count} missing genes)")
    
    logger.info(f"")
    
    # Write combined TSV file
    combined_report_file = os.path.join(report_directory, "all_samples_gene_validation_report.tsv")
    with open(combined_report_file, 'w') as f:
        f.write("\n".join(combined_tsv_lines))


def check_no_alignment_and_refs_order(args, gene_median_lengths, gene_synonyms=None):
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
    
    # Check if refs_order contains entries
    if args.refs_order:
        # If refs_order is specified, no_alignment must be False
        if args.no_alignment:
            logger.error(f"{"[ERROR]:":10} Cannot specify --refs_order when --no_alignment is True. Please remove --no_alignment or remove --refs_order.")
            utils.exit_program()
        
        # Get available order directories
        data_dir = os.path.join(os.path.dirname(__file__), 'data')
        order_genomes_dir = os.path.join(data_dir, 'order_genomes')
        
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
        order_refs = get_ref_gene_seqrecords_from_orders(args.refs_order, gene_median_lengths, gene_synonyms)
    
    # Check if custom reference folder is specified
    if args.custom_refs_folder:
        # If custom folder is specified, no_alignment must be False
        if args.no_alignment:
            logger.error(f"{"[ERROR]:":10} Cannot specify --custom_refs_folder when --no_alignment is True. Please remove --no_alignment or remove --custom_refs_folder.")
            utils.exit_program()
        
        try:
            logger.info(f"{"[INFO]:":10} Using custom reference genomes from: {args.custom_refs_folder}")
            custom_refs = get_ref_gene_seqrecords_from_custom_folder(args.custom_refs_folder, gene_median_lengths, gene_synonyms)
        except Exception as e:
            logger.error(f"{"[ERROR]:":10} Error processing custom reference folder: {e}")
            utils.exit_program()
    
    # If no specific references are specified, use default reference genomes
    if not args.refs_order and not args.refs_folder:
        logger.info(f"{"[INFO]:":10} Using default reference genomes from data/reference_genomes_default")
        default_refs = get_ref_gene_seqrecords_from_default(gene_median_lengths, gene_synonyms)
    
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


def get_ref_gene_seqrecords_from_orders(orders, gene_median_lengths, gene_synonyms=None):
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

    # Get path to the order genomes data directory
    data_dir = os.path.join(os.path.dirname(__file__), 'data', 'order_genomes')
    
    # Track gene copies per file for ID generation
    gene_copy_counts = {}  # (mapped_gene_name, gbk_file) -> copy_count
    
    for order in orders:
        order_dir = os.path.join(data_dir, order)
        gbk_files = glob.glob(os.path.join(order_dir, '*.gb*'))
        
        logger.debug(f"{"[DEBUG]:":10} Found {len(gbk_files)} GenBank files in {order} directory")
        
        for gbk_file in gbk_files:
            try:
                # Parse GenBank file (handle both gzipped and uncompressed files)
                if gbk_file.endswith('.gz'):
                    with gzip.open(gbk_file, 'rt') as f:
                        records = list(SeqIO.parse(f, "genbank"))
                else:
                    records = list(SeqIO.parse(gbk_file, "genbank"))
                
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
    logger.debug(f"{"[DEBUG]:":10} Extracted {total_sequences} gene sequences for {total_genes} unique genes from order-specific GenBank files")
    
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
            # Parse GenBank file (handle both gzipped and uncompressed files)
            if gbk_file.endswith('.gz'):
                with gzip.open(gbk_file, 'rt') as f:
                    records = list(SeqIO.parse(f, "genbank"))
            else:
                records = list(SeqIO.parse(gbk_file, "genbank"))
            
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


def get_ref_gene_seqrecords_from_default(gene_median_lengths, gene_synonyms=None):
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

    # Get path to the default reference genomes directory
    data_dir = os.path.join(os.path.dirname(__file__), 'data', 'reference_genomes_default')
    
    # Track gene copies per file for ID generation
    gene_copy_counts = {}  # (mapped_gene_name, gbk_file) -> copy_count
    gbk_files = glob.glob(os.path.join(data_dir, '*.gb*'))
    
    logger.debug(f"{"[DEBUG]:":10} Found {len(gbk_files)} GenBank files in default reference data directory")
    
    for gbk_file in gbk_files:
        try:
            # Parse GenBank file (handle both gzipped and uncompressed files)
            if gbk_file.endswith('.gz'):
                with gzip.open(gbk_file, 'rt') as f:
                    records = list(SeqIO.parse(f, "genbank"))
            else:
                records = list(SeqIO.parse(gbk_file, "genbank"))
            
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
            
        gene_results = result.get('gene_results', {})
        gene_cds_info = result.get('gene_cds_info', {})
        
        for gene_name, gene_info in gene_results.items():
            # Only create alignments for genes that have CDS info and are in reference
            if gene_name in gene_cds_info:
                if gene_name in ref_gene_seqrecords:
                    cds_info_list = gene_cds_info[gene_name]
                    ref_sequences = ref_gene_seqrecords[gene_name]

                    alignment_tasks.append({
                        'sample_name': sample_name,
                        'gene_name': gene_name,
                        'cds_info_list': cds_info_list,
                        'ref_sequences': ref_sequences,
                        'outdir_alignments': sample_alignments_dir,  # Use sample-specific directory
                        'threads': threads
                    })
                else:
                    logger.warning(f"{"[WARNING]:":10} No reference sequences found for gene {gene_name}")
    
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
                raise result[0]

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
            
        gene_results = result.get('gene_results', {})
        gene_rRNA_info = result.get('gene_rRNA_info', {})
        
        for gene_name, gene_info in gene_results.items():

            # Only create alignments for genes that have rRNA info and are in reference
            if gene_name in gene_rRNA_info:

                if gene_name in ref_gene_seqrecords:
                    rRNA_info_list = gene_rRNA_info[gene_name]
                    ref_sequences = ref_gene_seqrecords[gene_name]

                    alignment_tasks.append({
                        'sample_name': sample_name,
                        'gene_name': gene_name,
                        'rRNA_info_list': rRNA_info_list,
                        'ref_sequences': ref_sequences,
                        'outdir_alignments': sample_alignments_dir,  # Use sample-specific directory
                        'threads': threads
                    })
                else:
                    logger.debug(f"{"[DEBUG]:":10} No reference sequences found for rRNA gene {gene_name}")
    
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
                raise result[0]

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
            
        gene_results = result.get('gene_results', {})
        gene_tRNA_info = result.get('gene_tRNA_info', {})
        
        for gene_name, gene_info in gene_results.items():

            # Only create alignments for genes that have tRNA info and are in reference
            if gene_name in gene_tRNA_info:

                if gene_name in ref_gene_seqrecords:


                    tRNA_info_list = gene_tRNA_info[gene_name]
                    ref_sequences = ref_gene_seqrecords[gene_name]

                    alignment_tasks.append({
                        'sample_name': sample_name,
                        'gene_name': gene_name,
                        'tRNA_info_list': tRNA_info_list,
                        'ref_sequences': ref_sequences,
                        'outdir_alignments': sample_alignments_dir,  # Use sample-specific directory
                        'threads': threads
                    })
                else:
                    logger.debug(f"{"[DEBUG]:":10} No reference sequences found for tRNA gene {gene_name}")
    
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
                raise result[0]

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
                try:
                    translated = utils.pad_seq(seqrecord)[0].translate()
                    if "*" in str(translated.seq)[:-1]:  # Exclude the last position
                        has_internal_stops = True
                        break
                except:
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
                    # subprocess.run(['trimal', '-backtrans', temp_nucleotide_path, '-in', aligned_protein_file, '-out', alignment_file, ], 
                    #              stderr=subprocess.PIPE, check=True)
                    
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
    gene_types = {}      # gene_name -> 'CDS' or 'rRNA'
    
    for sample_name, result in all_sample_results.items():
        if 'error' in result:
            continue
            
        gene_cds_info = result.get('gene_cds_info', {})
        gene_rRNA_info = result.get('gene_rRNA_info', {})
        gene_tRNA_info = result.get('gene_tRNA_info', {})
        
        # Process CDS genes
        for gene_name, cds_info_list in gene_cds_info.items():
            if gene_name not in gene_sequences:
                gene_sequences[gene_name] = []
                gene_types[gene_name] = 'CDS'
            
            for i, cds_info in enumerate(cds_info_list):
                seq = cds_info['cds_seq']
                seqrecord = SeqRecord(seq=seq, id=f"{sample_name}_{gene_name}_copy_{i+1}", description=f"CDS from {sample_name}")
                gene_sequences[gene_name].append(seqrecord)
        
        # Process rRNA genes
        for gene_name, rRNA_info_list in gene_rRNA_info.items():
            if gene_name not in gene_sequences:
                gene_sequences[gene_name] = []
                gene_types[gene_name] = 'rRNA'
            
            for i, rRNA_info in enumerate(rRNA_info_list):
                seq = rRNA_info['rRNA_seq']
                seqrecord = SeqRecord(seq=seq, id=f"{sample_name}_{gene_name}_copy_{i+1}", description=f"rRNA from {sample_name}")
                gene_sequences[gene_name].append(seqrecord)
        
        # Process tRNA genes
        for gene_name, tRNA_info_list in gene_tRNA_info.items():
            if gene_name not in gene_sequences:
                gene_sequences[gene_name] = []
                gene_types[gene_name] = 'tRNA'
            
            for i, tRNA_info in enumerate(tRNA_info_list):
                seq = tRNA_info['tRNA_seq']
                seqrecord = SeqRecord(seq=seq, id=f"{sample_name}_{gene_name}_copy_{i+1}", description=f"tRNA from {sample_name}")
                gene_sequences[gene_name].append(seqrecord)
    
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
                raise result[0]

    logger.debug(f"{"[DEBUG]:":10} Per-gene alignment generation complete")


def align_genes(all_sample_results, ref_gene_seqrecords, output_directory, pool_size, threads, refs_order):
    """
    Generate alignments for annotated genes with reference genes.
    
    Args:
        all_sample_results (dict): Results from gene checking (contains gene_results for each sample)
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


def linearize_genome_upstream_gene(gbk_file, fasta_file, output_dir, sample_name, linearize_gene='psbA'):
    """
    Linearize genome upstream of specified gene and write new fasta file.
    
    Args:
        gbk_file (str): Path to annotated GenBank file
        fasta_file (str): Path to original fasta file
        output_dir (str): Output directory for linearized fasta
        sample_name (str): Sample name for output file
        linearize_gene (str): Gene name to use for linearization (default: 'psbA')
        
    Returns:
        str: Path to the new linearized fasta file
    """
    try:
        # Read the annotated GenBank file
        gbk_io, adjusted_gbk_content = utils.build_adjusted_genbank_io(gbk_file, fasta_file, logger=logger)
        record = next(SeqIO.parse(gbk_io, 'genbank'))
        
        # Find the specified gene
        gene_start = None
        for feature in record.features:
            if feature.type == 'CDS' and 'gene' in feature.qualifiers:
                if feature.qualifiers['gene'][0] == linearize_gene:
                    gene_start = feature.location.start
                    break
        
        if gene_start is None:
            logger.warning(f"{"[WARNING]:":10} {linearize_gene} gene not found in {sample_name}, using original sequence")
            return fasta_file
        
        # Read original fasta sequence
        with open(fasta_file, 'r') as handle:
            fasta_record = next(SeqIO.parse(handle, 'fasta'))
        
        # Linearize upstream of the specified gene
        sequence = fasta_record.seq
        linearized_sequence = sequence[gene_start:] + sequence[:gene_start]
        
        # Create new fasta record
        linearized_record = SeqRecord(
            seq=linearized_sequence,
            id=fasta_record.id,
            description=f"Linearized upstream of {linearize_gene} (original position: {gene_start})"
        )
        
        # Write linearized fasta file
        linearized_fasta = os.path.join(output_dir, f"{sample_name}_linearized.fasta")
        with open(linearized_fasta, 'w') as handle:
            SeqIO.write(linearized_record, handle, 'fasta')
        
        logger.debug(f"{"[INFO]:":10} Linearized {sample_name} upstream of {linearize_gene} (position {gene_start})")
        return linearized_fasta
        
    except Exception as e:
        logger.error(f"{"[ERROR]:":10} Failed to linearize {sample_name}: {str(e)}")
        return fasta_file


def annotate_genomes(genome_fasta_dir, output_directory, chloe_project_dir=None, chloe_script_path=None, linearize_gene='psbA'):
    """
    Annotate genome fasta files using chloe annotate command.
    
    Args:
        genome_fasta_dir (str): Path to directory containing genome fasta files
        output_directory (str): Path to directory to write annotated genomes
        logger: Logger instance for logging messages
        
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
    utils.log_separator(logger)
    logger.info(f"{"[INFO]:":10} Found {len(fasta_files)} fasta files to annotate:")
    for fasta_file in fasta_files:
        logger.info(f"{" ":10} {os.path.basename(fasta_file)}")
    utils.log_separator(logger)
    time.sleep(0.1)  # As tqdm output is not written to the queue, we need to sleep to ensure the logger messages have time to be written 

    annotated_genomes = dict({})
    
    # Process each fasta file
    for fasta_file in tqdm(fasta_files, desc=f"{"[INFO]:":10} {"Annotating genomes":<20}", file=sys.stdout):
        
        # Get the basename and filename prefix
        basename = os.path.basename(fasta_file)
        filename_prefix = os.path.splitext(basename)[0]
        
        # Remove additional extensions if present (e.g., .fasta.gz -> .fasta)
        if filename_prefix.endswith('.fasta'):
            filename_prefix = filename_prefix[:-6]  # Remove .fasta
        elif filename_prefix.endswith('.fa'):
            filename_prefix = filename_prefix[:-3]  # Remove .fa
        elif filename_prefix.endswith('.fas'):
            filename_prefix = filename_prefix[:-4]  # Remove .fas

        filename_prefix_no_dots = filename_prefix.replace('.', '_')

        # Initialise dict for sample data
        annotated_genomes[filename_prefix_no_dots] = {}
        
        # Define output gbk file path
        output_sample_dir = os.path.join(annotated_genomes_dir, filename_prefix_no_dots)
        os.makedirs(output_sample_dir, exist_ok=True)

        # Define output file paths
        output_sample_gbk = os.path.join(output_sample_dir, f"{filename_prefix}.chloe.gbk")
        output_sample_gff = os.path.join(output_sample_dir, f"{filename_prefix}.chloe.gff")
        output_sample_gbk_original = os.path.join(output_sample_dir, f"{filename_prefix}.chloe.original.gbk")
        output_sample_gff_original = os.path.join(output_sample_dir, f"{filename_prefix}.chloe.original.gff")
        output_sample_gbk_linearized = os.path.join(output_sample_dir, f"{filename_prefix}_linearized.chloe.gbk")
        output_sample_gff_linearized = os.path.join(output_sample_dir, f"{filename_prefix}_linearized.chloe.gff")
        output_sample_fasta_linearized = os.path.join(output_sample_dir, f"{filename_prefix}_linearized.fasta")

        # Remove any existing non-renamed files
        if os.path.exists(output_sample_gbk):
            os.remove(output_sample_gbk)
        if os.path.exists(output_sample_gff):
            os.remove(output_sample_gff)

        logger.debug(f"{"[DEBUG]:":10} Annotating {basename} -> {os.path.basename(output_sample_dir)}")

        # Determine what needs to be done
        has_original_annotation = utils.file_exists_and_not_empty(output_sample_gbk_original) and utils.file_exists_and_not_empty(output_sample_gff_original)
        has_linearized_annotation_and_fasta = utils.file_exists_and_not_empty(output_sample_gbk_linearized) and utils.file_exists_and_not_empty(output_sample_gff_linearized) and utils.file_exists_and_not_empty(output_sample_fasta_linearized)
        
        if has_linearized_annotation_and_fasta:
            logger.debug(f"{"[INFO]:":10} {filename_prefix} already completely processed with linearized sequence")
        elif has_original_annotation:
            logger.debug(f"{"[INFO]:":10} {filename_prefix} needs linearization and re-annotation")
        else:
            logger.debug(f"{"[INFO]:":10} {filename_prefix} needs initial annotation")

        # Resolve chloe project and script paths
        if (chloe_project_dir is None) != (chloe_script_path is None):
            logger.error(f"{"[ERROR]:":10} Both --chloe_project_dir and --chloe_script must be provided together.")
            utils.exit_program()

        if chloe_project_dir and chloe_script_path:
            chloe_project = f'--project={chloe_project_dir}'
            chloe_script = chloe_script_path
            logger.info(f"{"[INFO]:":10} Using user-specified chloe project and script")
        else:
            try:
                conda_prefix = os.environ['CONDA_PREFIX']
                logger.debug(f"{"[DEBUG]:":10} The CONDA_PREFIX is: {conda_prefix}")
            except KeyError:
                logger.debug(f"{"[DEBUG]:":10} CONDA_PREFIX environment variable is not set.")
                conda_prefix = None

            if not conda_prefix:
                logger.error(f"{"[ERROR]:":10} Cannot locate chloe project/script without CONDA_PREFIX or user-specified paths")
                utils.exit_program()

            chloe_project = f'--project={conda_prefix}/bin/chloe'
            chloe_script = f'{conda_prefix}/bin/chloe/chloe.jl'

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
            output_sample_dir,
            fasta_file
        ]

        # Execute annotation workflow
        if not has_original_annotation:
            # Do initial annotation
            logger.debug(f"{"[INFO]:":10} Performing initial annotation for {filename_prefix}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # Rename original annotation files to preserve original annotation
            if os.path.exists(output_sample_gbk):
                os.rename(output_sample_gbk, output_sample_gbk_original)
            if os.path.exists(output_sample_gff):
                os.rename(output_sample_gff, output_sample_gff_original)

            has_original_annotation = True

        # Linearize and re-annotate (if we have original annotation but no linearized)
        if has_original_annotation and not has_linearized_annotation_and_fasta:
            logger.debug(f"{"[INFO]:":10} Linearizing {filename_prefix} upstream of {linearize_gene}...")
            
            # Linearize the genome upstream of the specified gene
            output_sample_fasta_linearized = linearize_genome_upstream_gene(
                output_sample_gbk_original, 
                fasta_file, 
                output_sample_dir, 
                filename_prefix_no_dots,
                linearize_gene
            )

            if output_sample_fasta_linearized != fasta_file:
                # Re-annotate with linearized sequence

                # Re-run chloe with linearized fasta
                cmd[9] = output_sample_fasta_linearized  # Replace the fasta file in the command

                logger.debug(f"{"[DEBUG]:":10} Re-annotating {filename_prefix} with linearized sequence: {cmd}")
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
   
                logger.debug(f"{"[DEBUG]:":10} Re-annotation complete for {filename_prefix}")

        # Store the linearized files
        annotated_genomes[filename_prefix_no_dots]['gbk'] = output_sample_gbk_linearized
        annotated_genomes[filename_prefix_no_dots]['gff'] = output_sample_gff_linearized
        annotated_genomes[filename_prefix_no_dots]['fasta'] = output_sample_fasta_linearized

    utils.log_separator(logger)
    time.sleep(0.1)

    return annotated_genomes


def convert_gbk_to_embl(annotated_genomes_dict, output_directory, metadata_dict=None):
    """
    Convert GenBank files to EMBL format with sample-specific metadata.
    
    Args:
        annotated_genomes_dict (dict): Dictionary of annotated genome information
        output_directory (str): Output directory for EMBL files
        metadata_dict (dict, optional): Dictionary mapping fasta_filename to metadata containing:
            - project_id: Project identifier
            - locus_tag: Locus tag prefix
            - genus_species: Genus and species name
            - authors: Author list
            - manuscript_title: Manuscript title
            - journal_name: Journal name
            - complete_or_partial: Genome completeness status
    
    Returns:
        None
    """

    # Create EMBL output directory
    embl_output_dir = os.path.join(output_directory, '02_embl_files')
    os.makedirs(embl_output_dir, exist_ok=True)
    
    # Get default metadata for samples without specific metadata
    default_metadata = get_default_metadata()
 
    # Track metadata usage
    samples_with_metadata = 0
    samples_without_metadata = 0
    
    for sample_name, genome_info in annotated_genomes_dict.items():
        if 'error' in genome_info:
            logger.warning(f"{"[WARNING]:":10} Skipping {sample_name} - has error: {genome_info['error']}")
            continue
            
        gbk_file = genome_info.get('gbk')
        if not gbk_file or not os.path.exists(gbk_file):
            logger.warning(f"{"[WARNING]:":10} GenBank file not found for {sample_name}")
            continue

        fasta_file = genome_info.get('fasta')
        if not fasta_file or not os.path.exists(fasta_file):
            logger.warning(f"{"[WARNING]:":10} Fasta file not found for {sample_name}")
            continue
        
        # Get fasta filename for metadata lookup
        fasta_filename = os.path.basename(fasta_file).split('_linearized')[0]
        
        # Get metadata for this sample
        if metadata_dict and fasta_filename in metadata_dict:
            sample_metadata = metadata_dict[fasta_filename]
            samples_with_metadata += 1
            logger.debug(f"{"[DEBUG]:":10} Using metadata for {sample_name} (fasta: {fasta_filename})")
        else:
            sample_metadata = default_metadata
            samples_without_metadata += 1
            if metadata_dict:
                logger.debug(f"{"[DEBUG]:":10} No metadata found for {sample_name} (fasta: {fasta_filename}), using defaults")
            else:
                logger.debug(f"{"[DEBUG]:":10} Using default metadata for {sample_name}")
        
        # Extract metadata values
        locus_tag_prefix = sample_metadata['locus_tag']
        project_id = sample_metadata['project_id']
        genus_species = sample_metadata['genus_species']
        linear_or_circular = sample_metadata['linear_or_circular']
            
        try:
            # Read GenBank file using utils.build_adjusted_genbank_io   
            gbk_io, adjusted_gbk_content = utils.build_adjusted_genbank_io(gbk_file, fasta_file, logger=logger)
            record = next(SeqIO.parse(gbk_io, 'genbank'))
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
            
            # Modify features for EMBL output
            locus_counter = 1
            
            # First pass: collect all gene names and repeat regions to assign locus tags
            gene_locus_map = {}  # gene_name -> locus_tag
            repeat_locus_map = {}  # repeat_id -> locus_tag

            for feature in record.features:
                if feature.type in ['gene', 'CDS', 'tRNA', 'intron']:
                    # Get gene name from qualifiers
                    gene_name = None
                    if 'gene' in feature.qualifiers:
                        gene_name = feature.qualifiers['gene'][0]
                    elif 'product' in feature.qualifiers:
                        gene_name = feature.qualifiers['product'][0]
                    
                    if gene_name and gene_name not in gene_locus_map:
                        gene_locus_map[gene_name] = f"{locus_tag_prefix}_LOCUS{locus_counter}"
                        locus_counter += 1
                
                elif feature.type == 'repeat_region':
                    # For repeat regions, assign locus tags sequentially
                    repeat_id = f"repeat_{len(repeat_locus_map) + 1}"
                    repeat_locus_map[repeat_id] = f"{locus_tag_prefix}_LOCUS{locus_counter}"
                    locus_counter += 1
            
            # Second pass: apply locus tags and CDS-specific modifications
            for feature in record.features:
                # Remove /ID= and /parent= and /name qualifiers if present (not used in EMBL format)
                if 'ID' in feature.qualifiers:
                    del feature.qualifiers['ID']
                if 'parent' in feature.qualifiers:
                    del feature.qualifiers['parent']
                if 'name' in feature.qualifiers:
                    del feature.qualifiers['name']
                
                if feature.type in ['gene', 'CDS', 'tRNA', 'intron']:
                    # Get gene name from qualifiers
                    gene_name = None
                    if 'gene' in feature.qualifiers:
                        gene_name = feature.qualifiers['gene'][0]
                    elif 'product' in feature.qualifiers:
                        gene_name = feature.qualifiers['product'][0]
                    
                    if gene_name and gene_name in gene_locus_map:
                        feature.qualifiers['locus_tag'] = [gene_locus_map[gene_name]]
                
                elif feature.type == 'repeat_region':
                    # For repeat regions, assign locus tags based on their order in the features list
                    repeat_count = 0
                    for f in record.features:
                        if f.type == 'repeat_region':
                            repeat_count += 1
                            if f == feature:  # Found the current feature
                                repeat_id = f"repeat_{repeat_count}"
                                feature.qualifiers['locus_tag'] = [repeat_locus_map[repeat_id]]
                                break
                
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
            embl_filename = f"{sample_name}.embl"
            embl_filepath = os.path.join(embl_output_dir, embl_filename)
            
            # Write EMBL file first pass
            with open(embl_filepath, 'w') as handle:
                SeqIO.write(record, handle, 'embl')

            # Convert EMBL file to ENA template format
            convert_embl_to_ena_template(embl_filepath, sample_metadata, record_seq_length, record.id)
            
            logger.debug(f"{"[DEBUG]:":10} Converted {sample_name} to EMBL format: {embl_filepath}")
            
        except Exception as e:
            logger.error(f"{"[ERROR]:":10} Failed to convert {sample_name} to EMBL format: {str(e)}\n{traceback.format_exc()}")
            continue
    
    logger.info(f"{"[INFO]:":10} EMBL conversion complete. Files written to: {embl_output_dir}")
    
    # Log metadata usage summary
    if metadata_dict:
        logger.debug(f"{"[DEBUG]:":10} Metadata usage: {samples_with_metadata} samples with metadata, {samples_without_metadata} samples with defaults")
    else:
        logger.debug(f"{"[DEBUG]:":10} All {samples_with_metadata + samples_without_metadata} samples used default metadata")
    
    utils.log_separator(logger)


def convert_embl_to_ena_template(embl_filepath, sample_metadata, record_seq_length, record_id):
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
    # TODO: Implement EMBL to ENA template conversion
    # This stub function should:
    # 1. Read the EMBL file using BioPython
    # 2. Extract sequence and annotation information
    # 3. Format the data according to ENA submission requirements
    # 4. Incorporate the provided sample metadata
    # 5. Write the ENA template file

  
    logger.debug(f"{"[DEBUG]:":10} Converting EMBL file to ENA template: {embl_filepath}")
    logger.debug(f"{"[DEBUG]:":10} Using metadata: {sample_metadata}")

    lines_to_write = []
    space_line_template = 'XX'

    id_line_template = f"ID   XXX; XXX; {sample_metadata['linear_or_circular']}; XXX; XXX; XXX; XXX."
    lines_to_write.append(id_line_template)
    lines_to_write.append(space_line_template)

    ac_line_template = f"AC   XXX;"
    lines_to_write.append(ac_line_template)
    lines_to_write.append(space_line_template)

    ac_line_template_2 = f"AC * _{record_id};"
    lines_to_write.append(ac_line_template_2)
    lines_to_write.append(space_line_template)

    pr_line_template = f"PR   Project:{sample_metadata['project_id']};"
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
    with open(embl_filepath.replace('.embl', '.ena.embl'), 'w') as handle:
        for line in lines_to_write:
            handle.write(line + '\n')   

    logger.debug(f"{"[DEBUG]:":10} ENA template file written to: {embl_filepath.replace('.embl', '.ena.embl')}")
    logger.debug(f"{"[DEBUG]:":10} ENA template file written to: {embl_filepath.replace('.embl', '.ena.embl')}")
   
    
def query_intergenic_regions(annotated_genomes_dict, output_directory, min_intergenic_length=50, blast_evalue=1e-10, 
                             debug_intergenic=False, max_blast_hits=1, pool_size=1, threads=1, log_queue=None):
    """
    Query intergenic regions from annotated genomes.
    
    Args:
        annotated_genomes_dict (dict): Dictionary of annotated genome information
        output_directory (str): Output directory for results
        min_intergenic_length (int): Minimum length of intergenic region to analyze
        blast_evalue (float): BLAST E-value threshold
    
    Returns:
        None
    """
    logger.info(f"{"[INFO]:":10} Querying intergenic regions from annotated genomes...")

    # Create output directory for intergenic analysis
    intergenic_output_dir = os.path.join(output_directory, "04_intergenic_analysis")
    os.makedirs(intergenic_output_dir, exist_ok=True)

    # Path to BLAST database
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
            
        sample_data = {
            'sample_name': sample_name,
            'genome_info': genome_info,
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
                write_blast_report(existing_results, combined_report_file, "ALL_GENOMES")
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
        write_blast_report(all_results, combined_report_file, "ALL_GENOMES")
        logger.info(f"{"[INFO]:":10} Combined report written to: {combined_report_file}")
    else:
        logger.warning(f"{"[WARNING]:":10} No BLAST results found for any genome")

    logger.info(f"{"[INFO]:":10} Intergenic region analysis complete. Results in: {intergenic_output_dir}")
    utils.log_separator(logger)


def extract_intergenic_regions(gbk_file, fasta_file, min_length=50, debug_intergenic=False, logger=None):
    """
    Extract intergenic regions from a GenBank file.
    
    Args:
        gbk_file (str): Path to GenBank file
        fasta_file (str): Path to corresponding FASTA file
        min_length (int): Minimum length of intergenic region to extract
    
    Returns:
        list: List of dictionaries containing intergenic region information
    """
    intergenic_regions = []
    
    try:
        # Parse GenBank file
        gbk_io, adjusted_gbk_content = utils.build_adjusted_genbank_io(
            gbk_file_path=gbk_file,
            fasta_file_path=fasta_file,
            logger=logger,
        )

        record = SeqIO.read(gbk_io, "genbank")
        
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
                    
    except Exception as e:
        if logger:
            logger.error(f"{"[ERROR]:":10} Error extracting intergenic regions from {gbk_file}: {str(e)}")
            logger.error(traceback.format_exc())
    
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
        genome_info = sample_data['genome_info']
        blast_db_path = sample_data['blast_db_path']
        min_intergenic_length = sample_data['min_intergenic_length']
        blast_evalue = sample_data['blast_evalue']
        debug_intergenic = sample_data['debug_intergenic']
        max_blast_hits = sample_data['max_blast_hits']
        threads = sample_data['threads']
        intergenic_output_dir = sample_data['intergenic_output_dir']
        
        # Get GenBank file path
        gbk_file = genome_info.get('gbk')
        if not gbk_file or not os.path.exists(gbk_file):
            return False, f"GenBank file not found for {sample_name}"

        # Get FASTA file path
        fasta_file = genome_info.get('fasta')
        if not fasta_file or not os.path.exists(fasta_file):
            return False, f"FASTA file not found for {sample_name}"

        # Extract intergenic regions
        intergenic_regions = extract_intergenic_regions(gbk_file, fasta_file, min_intergenic_length, debug_intergenic, worker_logger)
        
        if not intergenic_regions:
            return False, (ValueError(f"No intergenic regions found for {sample_name}"), traceback.format_exc())

        # BLAST intergenic regions
        blast_results = blast_intergenic_regions(intergenic_regions, blast_db_path, blast_evalue, sample_name, max_blast_hits, threads, worker_logger)
        
        # Write individual genome report
        genome_report_file = os.path.join(intergenic_output_dir, f"{sample_name}_intergenic_blast_results.tsv")
        write_blast_report(blast_results, genome_report_file, sample_name, worker_logger)
        
        return True, (blast_results, genome_report_file)
        
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
        try:
            # Create temporary FASTA file for the intergenic region
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
                temp_fasta.write(f">{region['region_id']}\n{region['sequence']}\n")
                temp_fasta_path = temp_fasta.name
            
            # Create temporary output file for BLAST results
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_output:
                temp_output_path = temp_output.name
            
            # Run BLAST
            blast_cmd = [
                'blastn',
                '-db', blast_db_path,
                '-query', temp_fasta_path,
                '-out', temp_output_path,
                '-evalue', str(evalue_threshold),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle',
                '-max_target_seqs', str(max_blast_hits),
                '-num_threads', str(threads)
            ]
            
            result = subprocess.run(blast_cmd, capture_output=True, text=True, check=True)
            
            # Parse BLAST results
            with open(temp_output_path, 'r') as f:
                for line in f:
                    if line.strip():
                        fields = line.strip().split('\t')
                        if len(fields) >= 12:
                            sseqid = fields[1]
                            hit_origin = '_'.join(sseqid.split('[')[0].split('_')[:-1])
                            hit_origin = '_'.join(hit_origin.split('_')[:-1])
                            hit_gene = sseqid.split('[')[-1].rstrip(']')
                            hit_gene = '_'.join(hit_gene.split('_')[1:])
                                                
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
                            blast_results.append(blast_result)
            
            # Clean up temporary files
            os.unlink(temp_fasta_path)
            os.unlink(temp_output_path)
            
        except subprocess.CalledProcessError as e:
            if logger:
                logger.error(f"{"[ERROR]:":10} BLAST failed for region {region['region_id']}: {str(e)}")
                logger.error(f"{"[ERROR]:":10} BLAST stderr: {e.stderr}")
        except Exception as e:
            if logger:
                logger.error(f"{"[ERROR]:":10} Error processing region {region['region_id']}: {str(e)}")
    
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


def write_blast_report(blast_results, output_file, sample_name, logger=None):
    """
    Write BLAST results to TSV file.
    
    Args:
        blast_results (list): List of BLAST result dictionaries
        output_file (str): Output file path
        sample_name (str): Name of the sample or "ALL_GENOMES" for combined report
    """
    if not blast_results:
        if logger:
            logger.warning(f"{"[WARNING]:":10} No BLAST results to write for {sample_name}")
        return
    
    try:
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
        
        if logger:
            logger.debug(f"{"[DEBUG]:":10} BLAST report written for {sample_name}: {output_file}")
        
    except Exception as e:
        if logger:
            logger.error(f"{"[ERROR]:":10} Error writing BLAST report for {sample_name}: {str(e)}")
            logger.error(traceback.format_exc())


def parse_metadata_tsv(metadata_file_path, logger):
    """
    Parse and validate the metadata TSV file for EMBL conversion.
    
    This function reads a TSV file containing sample metadata and validates:
    1. File format and required columns
    2. Data completeness and format
    3. Consistency of the data
    
    Args:
        metadata_file_path (str): Path to the metadata TSV file
        logger: Logger instance for logging messages and errors
        
    Returns:
        dict: Dictionary mapping fasta_filename to metadata dictionary containing:
            - project_id: Project identifier
            - locus_tag: Locus tag prefix
            - genus_species: Genus and species name
            - linear_or_circular: Whether the genome is linear or circular
            
    Raises:
        SystemExit: If the file cannot be parsed or contains invalid data
    """
    required_columns = [
        'fasta_filename', 'project_id', 'locus_tag', 'genus_species', 'linear_or_circular'
    ]
    
    try:
        # Read the TSV file
        metadata_df = pd.read_csv(metadata_file_path, sep='\t', dtype=str)
        
        # Check for required columns
        missing_columns = set(required_columns) - set(metadata_df.columns)
        if missing_columns:
            logger.error(f"{"[ERROR]:":10} Missing required columns in metadata TSV: {missing_columns}")
            logger.error(f"{"":10} Required columns: {required_columns}")
            utils.exit_program()
        
        # Check for empty dataframe
        if metadata_df.empty:
            logger.error(f"{"[ERROR]:":10} Metadata TSV file is empty")
            utils.exit_program()
        
        # Check for duplicate fasta filenames
        duplicates = metadata_df[metadata_df.duplicated(['fasta_filename'], keep=False)]
        if not duplicates.empty:
            logger.error(f"{"[ERROR]:":10} Duplicate fasta_filename entries found:")
            for _, row in duplicates.iterrows():
                logger.error(f"{"":10} {row['fasta_filename']}")
            utils.exit_program()
        
        # Check for missing values in required fields
        required_fields = ['fasta_filename', 'project_id', 'locus_tag', 'genus_species', 'linear_or_circular']
        for field in required_fields:
            missing_values = metadata_df[metadata_df[field].isna() | (metadata_df[field] == '')]
            if not missing_values.empty:
                logger.error(f"{"[ERROR]:":10} Missing values in required field '{field}':")
                for _, row in missing_values.iterrows():
                    logger.error(f"{"":10} {row['fasta_filename']}")
                utils.exit_program()
        
        # Convert to dictionary
        metadata_dict = {}
        for _, row in metadata_df.iterrows():
            fasta_filename = row['fasta_filename'].strip()
            metadata_dict[fasta_filename] = {
                'project_id': row['project_id'].strip(),
                'locus_tag': row['locus_tag'].strip(),
                'genus_species': row['genus_species'].strip(),
                'linear_or_circular': row['linear_or_circular'].strip()
            }
        
        logger.info(f"{"[INFO]:":10} Successfully parsed metadata TSV file: {len(metadata_dict)} samples")
        logger.debug(f"{"[DEBUG]:":10} Sample entries: {list(metadata_dict.keys())[:5]}")
        
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


def get_default_metadata():
    """
    Get default metadata values for samples without metadata.
    
    Returns:
        dict: Dictionary containing default metadata values
    """
    return {
        'project_id': 'UNKNOWN_PROJECT',
        'locus_tag': 'DEFAULT_TAG',
        'genus_species': 'Unknown species',
        'linear_or_circular': 'linear_or_circular'
    }


def check_metadata_coverage(annotated_genomes_dict, metadata_dict, logger):
    """
    Check metadata coverage for fasta files and provide information messages.
    
    Args:
        annotated_genomes_dict (dict): Dictionary of annotated genome information
        metadata_dict (dict): Dictionary mapping fasta_filename to metadata
        logger: Logger instance for logging messages
        
    Returns:
        None
    """

    if not metadata_dict:
        return
    
    total_samples = len(annotated_genomes_dict)
    samples_with_metadata = 0
    samples_without_metadata = []
    
    for sample_name, genome_info in annotated_genomes_dict.items():
        if 'error' in genome_info:
            continue
            
        fasta_file = genome_info.get('fasta')
        if not fasta_file or not os.path.exists(fasta_file):
            continue
        
        fasta_filename = os.path.basename(fasta_file).split('_linearized')[0]
        if fasta_filename in metadata_dict:
            samples_with_metadata += 1
        else:
            samples_without_metadata.append(fasta_filename)
    
    logger.info(f"{"[INFO]:":10} Metadata coverage: {samples_with_metadata}/{total_samples} samples have metadata")
    
    if samples_without_metadata:
        logger.info(f"{"[INFO]:":10} Samples without metadata (will use defaults):\n")
        for fasta_filename in samples_without_metadata[:10]:  # Show first 10
            logger.info(f"{"":15} {fasta_filename}")
        if len(samples_without_metadata) > 10:
            logger.info(f"{"":15} ... and {len(samples_without_metadata) - 10} more")
        logger.info("")


def main(args):
    """Add sequences to an existing target file.
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

        # Set up global logger
        logger, log_queue, log_listener = utils.log_manager.setup(
            __name__, 'annotate_and_check', log_directory=args.log_directory
        )

        # Print arguments to screen and log:
        utils.print_arguments(args, logger, __version__)

        # Check for external dependencies:
        utils.check_dependencies(logger) 

        # Load gene median lengths from package resources
        gene_median_lengths = load_gene_median_lengths()

        # Load gene synonyms from package resources
        gene_synonyms = load_gene_synonyms()

        # Check no_alignment and refs_order
        ref_gene_seqrecords = check_no_alignment_and_refs_order(args, gene_median_lengths, gene_synonyms)
 
        # Annotate the genomes using chloë, honoring optional user-specified chloe paths
        annotated_genomes_dict = annotate_genomes(
            args.genome_fasta_dir,
            args.output_directory,
            getattr(args, 'chloe_project_dir', None),
            getattr(args, 'chloe_script', None),
            getattr(args, 'linearize_gene', 'psbA')
        )

        # Check genes and write reports
        all_sample_results = check_genes(gene_median_lengths, annotated_genomes_dict, args.min_length_percentage, args.max_length_percentage, 
                                         args.report_directory, log_queue, args.pool, gene_synonyms)

        # Parse metadata TSV file if provided
        metadata_dict = None
        logger.info(f"{"[INFO]:":10} Converting GenBank files to EMBL format...")
        if args.metadata_tsv:
            logger.info(f"{"[INFO]:":10} Parsing metadata TSV file: {args.metadata_tsv}")
            metadata_dict = parse_metadata_tsv(args.metadata_tsv, logger)
        else:
            logger.info(f"{"[INFO]:":10} No metadata TSV file provided, will use default values for EMBL conversion")
 
        # Check metadata coverage for fasta files
        check_metadata_coverage(annotated_genomes_dict, metadata_dict, logger)
 
        # Convert assembly gbk files to embl format
        convert_gbk_to_embl(annotated_genomes_dict, args.output_directory, metadata_dict=metadata_dict)

        # Generate alignments if not disabled
        if not args.no_alignment:
            align_genes(all_sample_results, ref_gene_seqrecords, args.output_directory, args.pool, args.threads, args.refs_order)

        # Query intergenic regions
        if not args.skip_intergenic_analysis:
            query_intergenic_regions(annotated_genomes_dict, args.output_directory, args.min_intergenic_length, args.blast_evalue, args.debug_intergenic, args.max_blast_hits, args.pool, args.threads, log_queue)
        else:
            logger.info(f"{"[INFO]:":10} Skipping intergenic region analysis as requested") 
 
    except Exception as e:
        utils.log_manager.handle_error(e, "main()")

    finally:
        # Log total completion time before cleaning up the logger
        utils.log_separator(logger)
        utils.log_completion_time(start_time, logger if ('logger' in globals() and logger) else None, label="PAV subcommand `annotate_and_check` completed")

        utils.log_manager.cleanup()