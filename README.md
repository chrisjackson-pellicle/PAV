# Plastid Annotation Validator (PAV)

A tool for annotating and validating angiosperm plastid genome annotations, with support for reference-based alignment, quality assessment, and characterisation of intergenic regions. 

## Overview

PAV is designed to process plastid genome assemblies, perform automated annotation using [Chloë](https://github.com/ian-small/Chloe.jl), validate gene annotations against reference sequences, and generate comprehensive reports. The tool includes features for reference-based alignment, and detailed quality assessment.

## Features

- **Automated Annotation**: Uses Chloë for plastid genome annotation
- **Genome Linearization**: Automatically linearizes genomes upstream of `psbA` gene
- **Reference-Based Validation**: Compares annotations against reference genomes from specified orders or default genome set
- **Quality Assessment**: Validates gene lengths, identifies any nternal stop codons, and checks for canonical start and stop codons
- **Alignment Generation**: Creates nucleotide (rRNA, tRNA) and protein alignments with reference sequences for CDS, rRNA, and tRNA genes
- **Comprehensive Reporting**: Generates detailed reports and statistics
- **EMBL and ENA Template Conversion**: Converts annotated GenBank records to EMBL and produces ENA submission-ready templates

## Installation

### Prerequisites

- Python 3.7+
- Chloë annotation tool
- MAFFT (for alignments)
- TrimAl (for backtranslation)
- BLAST+ (for querying intergenic regions)

### Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/chrisjackson-pellicle/PAV.git
   cd PAV
   ```

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Install Chloë** (if not already installed):
   ```bash
   # Follow Chloë installation instructions
   # If using a conda environment, ensure chloe and chloe references are in `${CONDA_PREFIX}/bin` or supply the chloe project dir and path to chloe.jl to the `pav` subcommand
   ```

4. **Install external tools**:
   ```bash
   # Install MAFFT
   conda install mafft

   # Install TrimAl
   conda install trimal

   # Install BLAST+
   conda install -c bioconda blast
   ```

## Usage

### Basic Usage

```bash
pav annotate_and_check \
  --input_dir /path/to/fasta/files \
  --output_dir /path/to/output \
```

### Command Line Options

#### Required Arguments
- `--input_dir`: Directory containing input FASTA files
- `--output_dir`: Directory for output files

#### Optional Arguments
- `--refs_order`: Reference order(s) to use for validation (can specify multiple)
- `--no_alignment`: Skip alignment generation
- `--min_length_percentage`: Minimum gene length percentage (default: 50)
- `--max_length_percentage`: Maximum gene length percentage (default: 200)
- `--pool`: Number of processes for multiprocessing (default: 1)
- `--threads`: Number of threads per process (default: 1)
 - `--metadata_tsv`: TSV file providing sample metadata for EMBL/ENA conversion
 - `--skip_intergenic_analysis`: Skip intergenic BLAST analysis
 - `--min_intergenic_length`: Minimum intergenic length to analyze (default: 0)
 - `--blast_evalue`: BLAST E-value threshold (default: 1e-10)
 - `--max_blast_hits`: Max BLAST hits to retain per region (default: 1)
 - `--debug_intergenic`: Write intergenic regions to FASTA for debugging

### Example Commands

#### Basic annotation and validation:
```bash
pav annotate_and_check \
  --input_dir genomes/ \
  --output_dir results/ \
```

#### Multiple reference orders:
```bash
pav annotate_and_check \
  --input_dir genomes/ \
  --output_dir results/ \
  --refs_order Alismatales Poales Arecales
```

#### Skip alignments:
```bash
pav annotate_and_check \
  --input_dir genomes/ \
  --output_dir results/ \
  --no_alignment
```

#### Custom parameters:
```bash
pav annotate_and_check \
  --input_dir genomes/ \
  --output_dir results/ \
  --refs_order Alismatales \
  --min_length_percentage 0.9 \
  --max_length_percentage 1.1 \
  --pool 4 \
  --threads 2
```

## Output Structure

The tool generates a structured output directory:

```
output_dir/
├── 00_logs_and_reports/
│   ├── logs/
│   │   └── annotate_and_check_<timestamp>.log          # Runtime log
│   └── reports/
│       ├── all_samples_gene_validation_report.tsv      # Combined gene length/translation checks
│       └── <sample_name>_gene_validation_report.tsv    # Per-sample reports
├── 01_annotated_genomes/
│   └── sample_name/
│       ├── <prefix>.chloe.original.gbk          # Original annotation (preserved)
│       ├── <prefix>.chloe.original.gff          # Original GFF (preserved)
│       ├── <prefix>_linearized.chloe.gbk        # Re-annotated after linearization
│       ├── <prefix>_linearized.chloe.gff        # Re-annotated after linearization
│       └── <prefix>_linearized.fasta            # Linearized sequence
│       └── <prefix>_linearized.chloe_intergenic_debug.fasta   # Optional (when --debug_intergenic)
├── 02_embl_files/
│   ├── <sample_name>.embl                        # EMBL format
│   └── <sample_name>.ena.embl                    # ENA template (derived from EMBL)
├── 03_alignments_with_refs/
│   ├── 01_per_sample_alignments/
│   │   └── <sample_name>/
│   │       ├── <sample>_<gene>_alignment.fasta          # Per-sample CDS alignment
│   │       ├── <sample>_<gene>_rRNA_alignment.fasta     # Per-sample rRNA alignment
│   │       └── <sample>_<gene>_tRNA_alignment.fasta     # Per-sample tRNA alignment
│   └── 02_per_gene_alignments/
│       ├── <gene>_all_samples_alignment.fasta           # All samples combined (nucleotide)
├── 04_intergenic_analysis/
    ├── <sample_name>_intergenic_blast_results.tsv       # Per-sample BLAST results
    └── combined_intergenic_blast_results.tsv            # Combined BLAST results

```

## Workflow

### 1. Genome Annotation
- Processes input FASTA files using Chloë
- Performs initial annotation on original sequences
- Linearizes genomes upstream of psbA gene
- Re-annotates linearized sequences

### 2. Annotation Validation
- Validates gene annotations against reference data
- Checks for gene length outliers, internal stop codons, and non-canonical start and stop codons
- Generates detailed reports

### 3. Alignment Generation
- Creates nucleotide alignments with reference sequences for rRNA and tRNA genes
- Generates protein alignments with backtranslated nucleotide alignments for CDS genes

### 4. EMBL and ENA Template Generation
- Adds locus tags and standardizes features for EMBL
- Builds EMBL template using metadata TSV (see below), ready for submission to ENA

## Metadata TSV Format

Provide a tab-separated file with the following required columns:

- `fasta_filename`
- `project_id`
- `locus_tag`
- `genus_species`
- `linear_or_circular` (must be either `linear` or `circular`)

Example (`metadata.tsv`):

```
fasta_filename	project_id	locus_tag	genus_species	linear_or_circular
sample1.fasta	PRJEB12345	ABC	Arabidopsis thaliana	circular
sample2.fasta	PRJEB98765	XYZ	Oryza sativa	linear
```

Notes:
- Only the columns above are required; defaults are used for anything else
- The topology is used in the ENA `ID` line and must match one of: `linear`, `circular`

## Gene Types Supported

PAV processes and validates three main types of plastid genes:

- **CDS (Coding Sequences)**: Protein-coding genes with translation validation
- **rRNA (Ribosomal RNA)**: Ribosomal RNA genes with nucleotide alignment
- **tRNA (Transfer RNA)**: Transfer RNA genes with nucleotide alignment

All gene types are extracted from reference GenBank files and used for validation and alignment generation.

## File Formats

### Input
- **FASTA**: Genome assemblies (.fasta, .fa, .fas)
- **Compressed**: Supports .gz compression

### Output
- **GenBank**: Annotated genomes in GenBank format
- **GFF**: Gene feature format files
- **EMBL**: European Molecular Biology Laboratory format
- **FASTA**: Aligned sequences and linearized genomes

## Reference Data

PAV uses reference genomes from the `data/reference_genomes_default` directory by default, or he  `data/order_genomes/` directory if orders are specified via the `--refs_order` flag. Each order subfolder contains:
- Reference genomes in GenBank format
- Taxonomic information for sequence identification

## Configuration


### Median Lengths
Reference gene lengths are stored in `data/plDNA_genes_median_lengths.csv` for validation.

## Troubleshooting

### Common Issues

1. **Chloë not found**: Ensure Chloë is properly installed and accessible
2. **MAFFT/TrimAl errors**: Check that external tools are installed and in PATH
3. **BLAST+ not found**: Install BLAST+ (`blastn`, `makeblastdb`) and ensure they are in your $PATH
4. **Metadata TSV columns**: Ensure the TSV has exactly the required columns listed above
5. **Memory issues**: Reduce `--pool` and `--threads` parameters


## Support

For issues and questions:
- Check the troubleshooting section
- Review the documentation
- Open an issue on GitHub

## Changelog

### Version 0.0.1

- Initial release
- Genome annotation with Chloë
- Reference-based validation
- Alignment generation
- Quality assessment and reporting
