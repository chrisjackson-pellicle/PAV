#!/usr/bin/env python3
"""
Script to extract genes and features from GenBank files in order_genomes directory.
Excludes inverted repeats, LSC, SSC, repeat_region, misc_feature, intron, gene, and source regions and writes features to FASTA format with proper naming.
"""

import os
import gzip
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import argparse


def parse_taxonomy_file(taxonomy_file):
    """Parse taxonomy file to get accession to species mapping."""
    accession_to_info = {}
    try:
        with open(taxonomy_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:  # Skip header if present
                for line in lines[1:]:
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        accession = parts[0]
                        species = parts[1]
                        family = parts[3]
                        order = parts[4]
                        accession_to_info[accession] = {
                            'species': species,
                            'family': family,
                            'order': order
                        }
    except Exception as e:
        print(f"Warning: Could not parse taxonomy file {taxonomy_file}: {e}")
    return accession_to_info


def extract_sequence_from_location(record, location):
    """Extract sequence from a feature location using BioPython's extract method."""
    try:
        # Use BioPython's extract method which handles complex locations correctly
        extracted_seq = location.extract(record.seq)
        return str(extracted_seq)
    except Exception as e:
        raise e
        # Fallback to manual extraction if extract fails
        if hasattr(location, 'parts'):  # CompoundLocation
            parts = []
            for part in location.parts:
                parts.append(str(record.seq[part.start:part.end]))
            return ''.join(parts)
        else:  # Simple location
            return str(record.seq[location.start:location.end])


def should_exclude_feature(feature):
    """Check if feature should be excluded (inverted repeats, source, LSC, SSC, repeat_region, misc_feature, intron, gene, etc.)."""
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


def process_genbank_file(gb_file, accession_to_info, order_name):
    """Process a single GenBank file and extract features."""
    features = []
    
    try:
        # Handle gzipped files
        if str(gb_file).endswith('.gz'):
            with gzip.open(gb_file, 'rt') as f:
                record = SeqIO.read(f, 'genbank')
        else:
            record = SeqIO.read(gb_file, 'genbank')
        
        # Get accession from filename
        accession = gb_file.stem.replace('.gb', '')
        
        # Get species info
        species_info = accession_to_info.get(accession, {})
        species = species_info.get('species', 'Unknown_species')
        family = species_info.get('family', 'Unknown_family')
        
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
                print(f"Warning: Could not extract sequence for feature in {gb_file}: {e}")
                continue
            
            # Create feature identifier
            feature_id = f"{order_name}_{family}_{species_clean}_{accession}"
            
            # Add feature type and gene name if available
            feature_desc = feature.type
            if 'gene' in feature.qualifiers:
                feature_desc += f"_{feature.qualifiers['gene'][0]}"
            elif 'product' in feature.qualifiers:
                feature_desc += f"_{feature.qualifiers['product'][0]}"
            
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
            
            features.append(seq_record)
            
    except Exception as e:
        print(f"Error processing {gb_file}: {e}")
    
    return features


def main():
    parser = argparse.ArgumentParser(description='Extract genes and features from GenBank files')
    parser.add_argument('--input_dir', default='plastid_annotation_validator/data/order_genomes',
                       help='Input directory containing order genomes')
    parser.add_argument('--output_file', default='extracted_features.fasta',
                       help='Output FASTA file')
    parser.add_argument('--min_length', type=int, default=10,
                       help='Minimum sequence length to include')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    output_file = args.output_file
    
    if not input_dir.exists():
        print(f"Error: Input directory {input_dir} does not exist")
        return
    
    all_features = []
    
    # Process each order directory
    for order_dir in sorted(input_dir.iterdir()):
        if not order_dir.is_dir() or order_dir.name.startswith('.'):
            continue
        
        order_name = order_dir.name
        print(f"Processing order: {order_name}")
        
        # Look for taxonomy file
        taxonomy_file = order_dir / 'plastid_genomes_taxonomy.tsv'
        accession_to_info = {}
        if taxonomy_file.exists():
            accession_to_info = parse_taxonomy_file(taxonomy_file)
        
        # Process GenBank files
        for gb_file in order_dir.glob('*.gb*'):
            if gb_file.name.endswith('.md') or gb_file.name.endswith('.tsv'):
                continue
            
            print(f"  Processing: {gb_file.name}")
            features = process_genbank_file(gb_file, accession_to_info, order_name)
            
            # Filter by minimum length
            features = [f for f in features if len(f.seq) >= args.min_length]
            all_features.extend(features)
    
    # Write to FASTA file
    print(f"Writing {len(all_features)} features to {output_file}")
    SeqIO.write(all_features, output_file, 'fasta')
    print(f"Done! Extracted {len(all_features)} features to {output_file}")


if __name__ == "__main__":
    main()
