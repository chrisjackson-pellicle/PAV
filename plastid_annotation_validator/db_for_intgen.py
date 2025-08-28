#!/usr/bin/env python3
"""
Script to extract genes and features from GenBank files in a single input directory.

This script processes GenBank files to extract gene sequences and features, excluding
inverted repeats, LSC, SSC, repeat_region, misc_feature, intron, gene, and source regions.
Features are written to FASTA format and then converted to a BLAST database.

The script expects GenBank files (.gb, .gb.gz, .gbk, .gbk.gz) in the input directory and creates
a BLAST database with sequences formatted as: Order_Family_genus_Species_accession_feature_type_[gene_name]

Taxonomy information is extracted directly from the GenBank file annotations.
"""

import os
import gzip
import re
import time
import traceback
import tempfile
import subprocess
import csv
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import glob

# Import PAV modules:
from plastid_annotation_validator.version import __version__
from plastid_annotation_validator import utils

# Initialise logger objects to None
logger = None
log_queue = None
log_listener = None

def extract_taxonomic_info(record):
    """Extract taxonomic information from GenBank record for ID generation.
    
    Args:
        record (SeqRecord): BioPython SeqRecord object containing GenBank annotations.
        
    Returns:
        dict: Dictionary with 'order', 'family', 'genus', 'species' keys.
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
    
    # Get order from taxonomy
    order = 'Unknown'
    # Look for order in taxonomy (usually ends with 'ales')
    for taxon in taxonomy:
        if taxon.endswith('ales'):
            order = taxon
            break
    
    return {
        'order': order,
        'family': family,
        'genus': genus,
        'species': species
    }


def extract_sequence_from_location(record, location):
    """Extract sequence from a feature location using BioPython's extract method.
    
    Args:
        record (SeqRecord): BioPython SeqRecord object containing the sequence.
        location (FeatureLocation): BioPython feature location object.
        
    Returns:
        str: Extracted DNA sequence as a string.
        
    Raises:
        Exception: If sequence extraction fails.
    """
    try:
        # Use BioPython's extract method which handles complex locations correctly
        extracted_seq = location.extract(record.seq)
        return str(extracted_seq)
    except Exception as e:
        raise e


def should_exclude_feature(feature):
    """Check if feature should be excluded based on type and qualifiers.
    
    Args:
        feature (SeqFeature): BioPython SeqFeature object to evaluate.
        
    Returns:
        bool: True if feature should be excluded, False otherwise.
    """
    # Exclude source features, LSC, SSC, repeat_region, misc_feature, intron, gene, and inverted repeats
    if feature.type in ['source', 'LSC', 'SSC', 'repeat_region', 'misc_feature', 'intron', 'gene']:
        return True
    
    # Check for inverted repeat in qualifiers
    if 'note' in feature.qualifiers:
        note = feature.qualifiers['note'][0].lower()
        if 'inverted repeat' in note or 'ir' in note:
            return True
    
    # Check for inverted repeat in gene name
    if 'gene' in feature.qualifiers:
        gene_name = feature.qualifiers['gene'][0].lower()
        if 'ir' in gene_name and ('inverted' in gene_name or 'repeat' in gene_name):
            return True
    
    return False


def process_genbank_file(gb_file):
    """Process a single GenBank file and extract features.
    
    Args:
        gb_file (str): Path to the GenBank file to process.
        
    Returns:
        list: List of SeqRecord objects representing extracted features.
    """
    features = []
    
    records = utils.parse_genbank_file(gb_file, logger)

    # Get accession from filename without any GenBank or gz extensions
    gb_path = Path(gb_file)
    accession = gb_path.stem.replace('.gb', '').replace('.gbk', '')

    for record in records:

        # Extract taxonomy information from GenBank record
        tax_info = extract_taxonomic_info(record)
        order = tax_info['order']
        family = tax_info['family']
        genus = tax_info['genus']
        species = tax_info['species']
    
        # Clean species name for FASTA header
        species_clean = species.replace(' ', '_').replace('var.', 'var').replace('subsp.', 'subsp')
    
        # Process each feature
        for feature in record.features:
            if should_exclude_feature(feature):
                continue
            
            # Skip features without location
            if not hasattr(feature, 'location') or feature.location is None:
                continue
            
            # Extract sequence
            try:
                sequence = extract_sequence_from_location(record, feature.location)
                if not sequence or len(sequence) == 0:
                    continue
            except Exception as e:
                if logger:
                    logger.warning(f"{"[WARNING]:":10} Could not extract sequence for feature in {gb_file}: {e}")
                continue
            
            # Create feature identifier
            feature_id = f"{order}_{family}_{genus}_{species_clean}_{accession}"

            # Add feature type and gene name if available
            feature_desc = ''
            if 'gene' in feature.qualifiers:
                feature_desc += f"{feature.qualifiers['gene'][0]}"
            elif 'product' in feature.qualifiers:
                feature_desc += f"{feature.qualifiers['product'][0]}"

            # Clean description for FASTA header
            feature_desc = re.sub(r'[^\w\-_]', '_', feature_desc)

            # Create FASTA header with feature type in square brackets
            fasta_header = f"{feature_id}_{feature.type}_[{feature_desc}]"

            # Create SeqRecord
            seq_record = SeqRecord(
                Seq(sequence),
                id=fasta_header,
                description=""
            )
            
            # Store source file information in the SeqRecord object
            seq_record.annotations['source_file'] = os.path.basename(gb_file)
            seq_record.annotations['accession'] = accession
            seq_record.annotations['feature_type'] = feature.type
            seq_record.annotations['gene_name'] = feature_desc if feature_desc else 'N/A'
            seq_record.annotations['order'] = order
            seq_record.annotations['family'] = family
            seq_record.annotations['genus'] = genus
            seq_record.annotations['species'] = species_clean

            features.append(seq_record)
        
    return features


def create_blast_database(fasta_file, output_dir, db_name="intgen_db"):
    """Create a BLAST database from a FASTA file.
    
    Args:
        fasta_file (str): Path to the input FASTA file.
        output_dir (str): Directory where the BLAST database will be created.
        db_name (str): Name for the BLAST database.
        
    Returns:
        str: Path to the created BLAST database.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Full path for the database
    db_path = os.path.join(output_dir, db_name)
    
    # Run makeblastdb command
    cmd = [
        'makeblastdb',
        '-in', fasta_file,
        '-dbtype', 'nucl',
        '-out', db_path,
        '-title', db_name,
        # '-parse_seqids'
    ]
    
    logger.debug(f"{"[DEBUG]:":10} Creating BLAST database: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        return db_path
    
    except subprocess.CalledProcessError as e:
        logger.error(f"{"[ERROR]:":10} Failed to create BLAST database: {e}")
        logger.error(f"{"[ERROR]:":10} STDOUT: {e.stdout}")
        logger.error(f"{"[ERROR]:":10} STDERR: {e.stderr}")
        raise
    except FileNotFoundError:
        logger.error(f"{"[ERROR]:":10} makeblastdb command not found. Please ensure BLAST+ is installed.")
        raise


def process_genbank_files(input_dir, min_length):
    """Process GenBank files in the input directory and extract features.
    
    Args:
        input_dir (str): Path to the input directory containing GenBank files.
        min_length (int): Minimum sequence length to include.
        
    Returns:
        list: List of SeqRecord objects representing extracted features.
    """

    # Ensure input directory exists
    if not os.path.exists(input_dir):
        logger.error(f"{"[ERROR]:":10} Input directory {input_dir} does not exist")
        utils.exit_program()
    
    all_features = []
    
    # Process GenBank files in the input directory
    gb_files = glob.glob(os.path.join(input_dir, "*.gb")) + \
                     glob.glob(os.path.join(input_dir, "*.gbk")) + \
                     glob.glob(os.path.join(input_dir, "*.gb.gz")) + \
                     glob.glob(os.path.join(input_dir, "*.gbk.gz"))
    
    if not gb_files:
        logger.error(f"{"[ERROR]:":10} No GenBank files found in {input_dir}")
        utils.exit_program()
    
    logger.info(f"{"[INFO]:":10} Found {len(gb_files)} GenBank files to process")
    
    for gb_file in gb_files:
        logger.info(f"{"[INFO]:":10} Processing: {os.path.basename(gb_file)}")
        features = process_genbank_file(gb_file)
        
        # Filter by minimum length
        features = [f for f in features if len(f.seq) >= min_length]
        all_features.extend(features)
        
    if not all_features:
        logger.error(f"{"[ERROR]:":10} No features extracted from GenBank files")
        utils.exit_program()
    
    utils.log_separator(logger)
    logger.info(f"{"[INFO]:":10} Extracted {len(all_features)} features from GenBank files")
    utils.log_separator(logger)
    
    return all_features


def write_blast_database(all_features, output_dir, db_name="intgen_db"):
    """Write features to a temporary FASTA file and create a BLAST database.
    
    Args:
        all_features (list): List of SeqRecord objects to write.
        output_dir (str): Directory where the BLAST database will be created.
        db_name (str): Name for the BLAST database.
        
    Returns:
        str: Path to the created BLAST database.
    """

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    # Create BLAST database output directory if it doesn't exist
    blast_db_dir = os.path.join(output_dir, f'01_{db_name}')
    os.makedirs(blast_db_dir, exist_ok=True)

    # Write to temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
        temp_fasta_path = temp_fasta.name
        logger.info(f"{"[INFO]:":10} Writing {len(all_features)} features to temporary FASTA file")
        SeqIO.write(all_features, temp_fasta_path, 'fasta')
        utils.log_separator(logger)
    
    try:
        # Create BLAST database
        db_path = create_blast_database(temp_fasta_path, blast_db_dir, db_name)
        logger.info(f"{"[INFO]:":10} Successfully created BLAST database at: {db_path}")
        logger.info(f"{"[INFO]:":10} Database contains {len(all_features)} sequences")
        utils.log_separator(logger)
        return db_path
    finally:
        # Clean up temporary FASTA file
        if os.path.exists(temp_fasta_path):
            os.unlink(temp_fasta_path)
            logger.debug(f"{"[DEBUG]:":10} Cleaned up temporary FASTA file: {temp_fasta_path}")


def write_tsv_report(all_features, report_directory):
    """Write a TSV report of extracted features.
    
    Args:
        all_features (list): List of SeqRecord objects to write.
        report_directory (str): Directory where the TSV report will be created.
    """
    report_path = os.path.join(report_directory, "extracted_features_report.tsv")
    
    with open(report_path, 'w', newline='') as f:
        fieldnames = [
            'Feature_ID', 'Sequence_Length', 'Feature_Type', 'Gene_Name', 'Product',
            'Order', 'Family', 'Genus', 'Species', 'Accession', 'Source_File'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        
        writer.writeheader()
        
        for feature in all_features:
            # Extract basic information
            seq_length = len(feature.seq)
            
            # Get information from stored annotations
            source_file = feature.annotations.get('source_file', 'unknown')
            accession = feature.annotations.get('accession', 'unknown')
            feature_type = feature.annotations.get('feature_type', 'unknown')
            gene_name = feature.annotations.get('gene_name', 'N/A')
            order = feature.annotations.get('order', 'unknown')
            family = feature.annotations.get('family', 'unknown')
            genus = feature.annotations.get('genus', 'unknown')
            species = feature.annotations.get('species', 'unknown')
            
            # Try to get product from description if available
            product = 'N/A'
            if hasattr(feature, 'description') and feature.description:
                product = feature.description
            
            writer.writerow({
                'Feature_ID': feature.id,
                'Sequence_Length': seq_length,
                'Feature_Type': feature_type,
                'Gene_Name': gene_name,
                'Product': product,
                'Order': order,
                'Family': family,
                'Genus': genus,
                'Species': species,
                'Accession': accession,
                'Source_File': source_file
            })
    
    logger.info(f"{"[INFO]:":10} TSV report written to: {report_path}")
    logger.info(f"{"[INFO]:":10} Report contains {len(all_features)} features")


def main(args):
    """Main function to process GenBank files and extract features.
    
    Args:
        args (argparse.Namespace): Parsed command line arguments containing input directory
            and other configuration options.
    """
    # Track wall-clock runtime for completion message
    start_time = time.time()
    
    try:
        global logger, log_queue, log_listener

        # Set up log and report directories
        utils.setup_log_and_report_directories(args)

        # Set up global logger
        logger, log_queue, log_listener = utils.log_manager.setup(
            __name__, 'db_for_intgen', log_directory=args.log_directory
        ) 
            
        # Print arguments to screen and log
        utils.print_arguments(args, logger, __version__)

        # Check for external dependencies:
        utils.check_dependencies(args, logger, entry='db_for_intgen')
        
        # Process GenBank files and extract features
        all_features = process_genbank_files(args.input_dir, args.min_length)
        
        # Create BLAST database
        write_blast_database(all_features, args.output_directory, "intgen_db")

        # Write TSV report
        write_tsv_report(all_features, args.report_directory)
        
    except Exception as e:
        utils.log_manager.handle_error(e, traceback.format_exc(), "db_for_intgen.main()")
    
    finally:
        # Log total completion time before cleaning up the logger
        utils.log_separator(logger)
        utils.log_completion_time(start_time, logger if ('logger' in globals() and logger) else None, 
                                  label="PAV subcommand `db_for_intgen` completed")
        utils.log_manager.cleanup()



