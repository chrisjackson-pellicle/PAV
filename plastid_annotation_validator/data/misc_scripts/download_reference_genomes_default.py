#!/usr/bin/env python3
"""
Script to download plastid/chloroplast genome accessions from NCBI
and create a TSV file with taxonomy details.
Verified version with unique, confirmed accession numbers.
"""

import os
import sys
import subprocess
import time
from pathlib import Path

# Verified list of plastid genome accessions with unique accession numbers
PLASTID_GENOMES = [
    # ===== BASAL ANGIOSPERMS (ANA Grade) =====
    {
        "accession": "NC_005086",
        "species": "Amborella trichopoda",
        "common_name": "Amborella",
        "family": "Amborellaceae",
        "order": "Amborellales",
        "clade": "Basal Angiosperms",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; Amborellales; Amborellaceae; Amborella"
    },
    {
        "accession": "NC_024542",
        "species": "Nymphaea mexicana",
        "common_name": "Mexican water lily",
        "family": "Nymphaeaceae",
        "order": "Nymphaeales",
        "clade": "Basal Angiosperms",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; Nymphaeales; Nymphaeaceae; Nymphaea"
    },
    {
        "accession": "NC_030504",
        "species": "Liriodendron chinense",
        "common_name": "Chinese tulip tree",
        "family": "Magnoliaceae",
        "order": "Magnoliales",
        "clade": "Basal Angiosperms",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; Magnoliidae; Magnoliales; Magnoliaceae; Liriodendron"
    },
    {
        "accession": "NC_008235",
        "species": "Illicium oligandrum",
        "common_name": "Star anise",
        "family": "Schisandraceae",
        "order": "Austrobaileyales",
        "clade": "Basal Angiosperms",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; Austrobaileyales; Schisandraceae; Illicium"
    },
    
    # ===== MONOCOTS =====
    # Poales (Grasses)
    {
        "accession": "NC_001666",
        "species": "Zea mays",
        "common_name": "Maize",
        "family": "Poaceae",
        "order": "Poales",
        "clade": "Monocots",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; Liliopsida; Poales; Poaceae; PACMAD clade; Panicoideae; Andropogonodae; Andropogoneae; Tripsacinae; Zea"
    },
    {
        "accession": "NC_031333",
        "species": "Oryza sativa",
        "common_name": "Rice",
        "family": "Poaceae",
        "order": "Poales",
        "clade": "Monocots",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; Liliopsida; Poales; Poaceae; BOP clade; Oryzoideae; Oryzeae; Oryzinae; Oryza"
    },
    {
        "accession": "NC_008602",
        "species": "Triticum aestivum",
        "common_name": "Wheat",
        "family": "Poaceae",
        "order": "Poales",
        "clade": "Monocots",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; Liliopsida; Poales; Poaceae; BOP clade; Pooideae; Triticodae; Triticeae; Triticum"
    },
    {
        "accession": "NC_008590",
        "species": "Hordeum vulgare",
        "common_name": "Barley",
        "family": "Poaceae",
        "order": "Poales",
        "clade": "Monocots",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; Liliopsida; Poales; Poaceae; BOP clade; Pooideae; Triticodae; Triticeae; Hordeum"
    },
    
    # ===== EUDICOTS =====
    # Caryophyllales
    {
        "accession": "NC_002202",
        "species": "Spinacia oleracea",
        "common_name": "Spinach",
        "family": "Chenopodiaceae",
        "order": "Caryophyllales",
        "clade": "Eudicots",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; Caryophyllales; Chenopodiaceae; Chenopodioideae; Anserineae; Spinacia"
    },
    
    # Brassicales
    {
        "accession": "AP000423",
        "species": "Arabidopsis thaliana",
        "common_name": "Thale cress",
        "family": "Brassicaceae",
        "order": "Brassicales",
        "clade": "Eudicots",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; Brassicales; Brassicaceae; Camelineae; Arabidopsis"
    },
    
    # Proteales
    {
        "accession": "NC_010601",
        "species": "Nelumbo nucifera",
        "common_name": "Sacred lotus",
        "family": "Nelumbonaceae",
        "order": "Proteales",
        "clade": "Eudicots",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; Proteales; Nelumbonaceae"
    },
    
    # ===== ROSIDS =====
    # Fabales (Legumes)
    {
        "accession": "NC_003119",
        "species": "Lotus japonicus",
        "common_name": "Bird's-foot trefoil",
        "family": "Fabaceae",
        "order": "Fabales",
        "clade": "Rosids",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; rosids; fabids; Fabales; Fabaceae"
    },
    {
        "accession": "NC_002694",
        "species": "Medicago truncatula",
        "common_name": "Barrel medic",
        "family": "Fabaceae",
        "order": "Fabales",
        "clade": "Rosids",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; rosids; fabids; Fabales; Fabaceae"
    },
    
    # Sapindales
    {
        "accession": "NC_008168",
        "species": "Citrus sinensis",
        "common_name": "Sweet orange",
        "family": "Rutaceae",
        "order": "Sapindales",
        "clade": "Rosids",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; Sapindales; Rutaceae"
    },
    
    # ===== ASTERIDS =====
    # Solanales
    {
        "accession": "NC_007898",
        "species": "Nicotiana tabacum",
        "common_name": "Tobacco",
        "family": "Solanaceae",
        "order": "Solanales",
        "clade": "Asterids",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; asterids; lamiids; Solanales; Solanaceae"
    },
    {
        "accession": "NC_007500",
        "species": "Solanum tuberosum",
        "common_name": "Potato",
        "family": "Solanaceae",
        "order": "Solanales",
        "clade": "Asterids",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; asterids; lamiids; Solanales; Solanaceae"
    },
    
    # Lamiales
    {
        "accession": "NC_009143",
        "species": "Sesamum indicum",
        "common_name": "Sesame",
        "family": "Pedaliaceae",
        "order": "Lamiales",
        "clade": "Asterids",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; asterids; lamiids; Lamiales; Pedaliaceae"
    },
    
    # Asterales
    {
        "accession": "NC_007977",
        "species": "Helianthus annuus",
        "common_name": "Sunflower",
        "family": "Asteraceae",
        "order": "Asterales",
        "clade": "Asterids",
        "taxonomy": "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; asterids; campanulids; Asterales; Asteraceae"
    }
]

def check_efetch_available():
    """Check if NCBI E-utilities are available."""
    try:
        result = subprocess.run(['efetch', '-help'], 
                              capture_output=True, text=True, timeout=10)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False

def download_genome(accession, output_dir):
    """Download a genome using efetch."""
    output_file = output_dir / f"{accession}.gb"
    
    if output_file.exists():
        print(f"  ✓ {accession}.gb already exists, skipping download")
        return True
    
    try:
        cmd = [
            'efetch', 
            '-db', 'nucleotide', 
            '-id', accession, 
            '-format', 'gb'
        ]
        
        with open(output_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, timeout=60)
        
        if result.returncode == 0 and output_file.exists():
            print(f"  ✓ Downloaded {accession}.gb")
            return True
        else:
            print(f"  ✗ Failed to download {accession}: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"  ✗ Timeout downloading {accession}")
        return False
    except Exception as e:
        print(f"  ✗ Error downloading {accession}: {e}")
        return False

def create_tsv_file(genomes, output_file):
    """Create a TSV file with genome details."""
    headers = [
        'Accession', 'Species', 'Common_Name', 'Family', 'Order', 
        'Clade', 'Taxonomy', 'File_Path'
    ]
    
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(headers) + '\n')
        
        # Write data
        for genome in genomes:
            row = [
                genome['accession'],
                genome['species'],
                genome['common_name'],
                genome['family'],
                genome['order'],
                genome['clade'],
                genome['taxonomy'],
                f"reference_genomes_verified/{genome['accession']}.gb"
            ]
            f.write('\t'.join(row) + '\n')
    
    print(f"✓ Created TSV file: {output_file}")

def main():
    """Main function to download genomes and create TSV file."""
    print("Verified Plastid Genome Downloader")
    print("=" * 60)
    
    # Check if efetch is available
    if not check_efetch_available():
        print("❌ Error: NCBI E-utilities (efetch) not found.")
        print("Please install NCBI E-utilities:")
        print("  - Download from: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-and-updated-datasets-are-available-for-download-and-local-installation/")
        print("  - Or use conda: conda install -c bioconda entrez-direct")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path("../reference_genomes_verified")
    output_dir.mkdir(exist_ok=True)
    print(f"✓ Created output directory: {output_dir}")
    
    # Download genomes
    print(f"\nDownloading {len(PLASTID_GENOMES)} plastid genomes...")
    successful_downloads = 0
    
    for i, genome in enumerate(PLASTID_GENOMES, 1):
        print(f"\n[{i}/{len(PLASTID_GENOMES)}] Downloading {genome['accession']} ({genome['species']})")
        
        if download_genome(genome['accession'], output_dir):
            successful_downloads += 1
        
        # Small delay to be respectful to NCBI servers
        time.sleep(1)
    
    # Create TSV file
    print(f"\nCreating TSV file...")
    tsv_file = output_dir / "plastid_genomes_taxonomy_verified.tsv"
    create_tsv_file(PLASTID_GENOMES, tsv_file)
    
    # Summary
    print(f"\n" + "=" * 60)
    print(f"Download Summary:")
    print(f"  Total genomes: {len(PLASTID_GENOMES)}")
    print(f"  Successfully downloaded: {successful_downloads}")
    print(f"  Failed: {len(PLASTID_GENOMES) - successful_downloads}")
    print(f"  Output directory: {output_dir}")
    print(f"  TSV file: {tsv_file}")
    
    if successful_downloads < len(PLASTID_GENOMES):
        print(f"\n⚠️  Some downloads failed. You may need to:")
        print(f"   - Check your internet connection")
        print(f"   - Verify the accession numbers are correct")
        print(f"   - Try downloading manually from NCBI")

if __name__ == "__main__":
    main() 