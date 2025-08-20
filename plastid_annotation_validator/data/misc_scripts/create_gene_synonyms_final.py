#!/usr/bin/env python3
"""
Script to create comprehensive gene_synonyms.txt with all gene names from GenBank files.
LOCATION: plastid_annotation_validator/data/misc_scripts/
UPDATED: Uses existing mappings from gene_synonyms.txt first, then fuzzy matching for new genes.
FIXED: Handle U to T variations for tRNA genes.
EXCLUDED: locus_tag values are not included in gene name extraction.
"""

import os
import gzip
import re
from difflib import SequenceMatcher
from collections import defaultdict, Counter
from Bio import SeqIO
import csv

def similarity(a, b):
    """Calculate similarity between two strings."""
    return SequenceMatcher(None, a.lower(), b.lower()).ratio()

def handle_u_to_t_variations(gene_name):
    """Generate U to T variations for tRNA genes."""
    variations = [gene_name]
    if 'trn' in gene_name.lower():
        # Convert U to T and vice versa
        if 'U' in gene_name:
            variations.append(gene_name.replace('U', 'T'))
        if 'T' in gene_name:
            variations.append(gene_name.replace('T', 'U'))
    return variations

def extract_gene_names_from_gb_file(file_path):
    """Extract all gene names from a GenBank file, excluding locus_tag."""
    gene_names = set()
    
    try:
        # Handle gzipped files
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as handle:
                records = SeqIO.parse(handle, "genbank")
                for record in records:
                    for feature in record.features:
                        if feature.type in ['gene', 'CDS', 'tRNA', 'rRNA']:
                            # Extract gene name from qualifiers, but exclude locus_tag
                            qualifiers_to_check = ['gene', 'product']
                            for qualifier in qualifiers_to_check:
                                if qualifier in feature.qualifiers:
                                    gene_name = feature.qualifiers[qualifier][0].strip()
                                    if gene_name and gene_name not in ['', 'hypothetical protein']:
                                        gene_names.add(gene_name)
        else:
            with open(file_path, 'r') as handle:
                records = SeqIO.parse(handle, "genbank")
                for record in records:
                    for feature in record.features:
                        if feature.type in ['gene', 'CDS', 'tRNA', 'rRNA']:
                            # Extract gene name from qualifiers, but exclude locus_tag
                            qualifiers_to_check = ['gene', 'product']
                            for qualifier in qualifiers_to_check:
                                if qualifier in feature.qualifiers:
                                    gene_name = feature.qualifiers[qualifier][0].strip()
                                    if gene_name and gene_name not in ['', 'hypothetical protein']:
                                        gene_names.add(gene_name)
                                        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
    
    return gene_names

def load_existing_mappings(synonyms_file):
    """Load existing gene mappings from gene_synonyms.txt file."""
    existing_mappings = {}
    if not os.path.exists(synonyms_file):
        return existing_mappings
    
    try:
        with open(synonyms_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue
                # Parse mapping: gene_name = standard_name
                if ' = ' in line:
                    parts = line.split(' = ', 1)
                    if len(parts) == 2:
                        gene_name = parts[0].strip()
                        standard_name = parts[1].strip()
                        existing_mappings[gene_name] = standard_name
    except Exception as e:
        print(f"Warning: Error loading existing mappings from {synonyms_file}: {e}")
    
    return existing_mappings

def load_standard_genes(csv_file):
    """Load standard gene names from CSV file."""
    standard_genes = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_name = row['Key'].strip()
            if gene_name:
                standard_genes.append(gene_name)
    return standard_genes

def find_best_match(gene_name, standard_genes, threshold=0.8):
    """Find the best matching standard gene name using fuzzy matching."""
    best_match = None
    best_score = 0
    
    # Check for exact matches (case insensitive)
    for std_gene in standard_genes:
        if gene_name.lower() == std_gene.lower():
            return std_gene, 1.0
    
    # Handle U to T variations for tRNA genes
    variations = handle_u_to_t_variations(gene_name)
    for variation in variations:
        for std_gene in standard_genes:
            if variation.lower() == std_gene.lower():
                return std_gene, 1.0
    
    # Fuzzy matching
    for std_gene in standard_genes:
        score = similarity(gene_name, std_gene)
        if score > best_score and score >= threshold:
            best_score = score
            best_match = std_gene
    
    return best_match, best_score

def main():
    # Paths (adjusted for script location in data/download_scripts/)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)  # Go up to data/ directory
    order_genomes_dir = os.path.join(base_dir, "order_genomes")
    reference_genomes_dir = os.path.join(base_dir, "reference_genomes_default")
    csv_file = os.path.join(base_dir, "plDNA_genes_median_lengths.csv")
    synonyms_file = os.path.join(base_dir, "gene_synonyms.txt")
    
    print("Loading existing gene mappings from gene_synonyms.txt...")
    existing_mappings = load_existing_mappings(synonyms_file)
    print(f"Loaded {len(existing_mappings)} existing gene mappings")
    
    print("Loading standard gene names from CSV...")
    standard_genes = load_standard_genes(csv_file)
    print(f"Loaded {len(standard_genes)} standard gene names")
    
    all_gene_names = set()
    file_count = 0
    
    print("Extracting gene names from reference genomes...")
    # Process reference genomes
    if os.path.exists(reference_genomes_dir):
        for filename in os.listdir(reference_genomes_dir):
            if filename.endswith('.gb'):
                file_path = os.path.join(reference_genomes_dir, filename)
                genes = extract_gene_names_from_gb_file(file_path)
                all_gene_names.update(genes)
                file_count += 1
                print(f"Processed {filename}: found {len(genes)} gene names")
    
    print("Extracting gene names from order genomes...")
    # Process order genomes
    if os.path.exists(order_genomes_dir):
        for order_name in os.listdir(order_genomes_dir):
            order_path = os.path.join(order_genomes_dir, order_name)
            if os.path.isdir(order_path):
                print(f"Processing order: {order_name}")
                for filename in os.listdir(order_path):
                    if filename.endswith('.gb.gz'):
                        file_path = os.path.join(order_path, filename)
                        genes = extract_gene_names_from_gb_file(file_path)
                        all_gene_names.update(genes)
                        file_count += 1
    
    print(f"\nExtracted {len(all_gene_names)} unique gene names from {file_count} files")
    
    # Create mappings
    mappings = {}
    exact_matches = 0
    fuzzy_matches = 0
    existing_used = 0
    no_matches = []
    
    print("Creating gene mappings...")
    for gene_name in sorted(all_gene_names):
        # Skip if it looks like a locus_tag (typically alphanumeric codes)
        if re.match(r'^[A-Z]+[0-9]+[A-Z]*[0-9]*$', gene_name):
            continue
        
        # First, check if we already have a mapping for this gene
        if gene_name in existing_mappings:
            mappings[gene_name] = existing_mappings[gene_name]
            existing_used += 1
        else:
            # If no existing mapping, try fuzzy matching
            best_match, score = find_best_match(gene_name, standard_genes)
            if best_match:
                mappings[gene_name] = best_match
                if score == 1.0:
                    exact_matches += 1
                else:
                    fuzzy_matches += 1
            else:
                no_matches.append(gene_name)
    
    # Create a backup of the existing file if it exists
    output_file = synonyms_file  # Use the dynamic path we already calculated
    backup_file = output_file + ".backup"
    
    if os.path.exists(output_file):
        print(f"Creating backup of existing file: {backup_file}")
        import shutil
        shutil.copy2(output_file, backup_file)
    
    # Group new mappings to add to the existing organized structure
    new_mappings = {}
    for gene_name, standard_name in mappings.items():
        if gene_name not in existing_mappings:
            new_mappings[gene_name] = standard_name
    
    # If we have existing mappings, preserve the organized structure
    if existing_mappings:
        print(f"Preserving existing organized structure with {len(existing_mappings)} mappings")
        print(f"Adding {len(new_mappings)} new mappings to the end")
        
        # Read the existing file and append new mappings
        with open(output_file, 'r') as f:
            existing_content = f.read()
        
        with open(output_file, 'w') as f:
            # Write existing content
            f.write(existing_content)
            
            # Add new mappings if any
            if new_mappings:
                f.write("\n# ================================================================================\n")
                f.write("# NEW MAPPINGS FOUND (require manual organization)\n")
                f.write("# ================================================================================\n\n")
                for gene_name in sorted(new_mappings.keys()):
                    f.write(f"{gene_name} = {new_mappings[gene_name]}\n")
            
            # Write unmatched genes as comments for manual review
            if no_matches:
                f.write("\n# ================================================================================\n")
                f.write("# UNMATCHED GENES REQUIRING MANUAL REVIEW\n")
                f.write("# ================================================================================\n")
                for gene_name in sorted(no_matches):
                    f.write(f"# {gene_name} = ?\n")
    else:
        # No existing mappings, create new file with basic structure
        with open(output_file, 'w') as f:
            f.write("# Gene synonyms mapping\n")
            f.write("# Format: extracted_gene_name = standard_gene_name\n")
            f.write(f"# Generated from {file_count} GenBank files\n")
            f.write(f"# Existing mappings: {existing_used}, Exact matches: {exact_matches}, Fuzzy matches: {fuzzy_matches}\n")
            f.write("# U to T variations handled for tRNA genes\n")
            f.write("# locus_tag values excluded\n")
            f.write("\n")
            
            # Write mappings
            for gene_name in sorted(mappings.keys()):
                f.write(f"{gene_name} = {mappings[gene_name]}\n")
            
            # Write unmatched genes as comments for manual review
            if no_matches:
                f.write("\n# Unmatched genes requiring manual review:\n")
                for gene_name in sorted(no_matches):
                    f.write(f"# {gene_name} = ?\n")
    
    print(f"\nResults:")
    print(f"- Total mapped genes: {len(mappings)}")
    print(f"  - Existing mappings used: {existing_used}")
    print(f"  - New exact matches: {exact_matches}")
    print(f"  - New fuzzy matches: {fuzzy_matches}")
    print(f"- New unmatched genes: {len(no_matches)}")
    print(f"- Output written to: {output_file}")
    if os.path.exists(backup_file):
        print(f"- Backup created: {backup_file}")

if __name__ == "__main__":
    main()
