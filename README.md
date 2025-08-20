# Plastid Annotation Validator (PAV)

A tool for annotating and validating angiosperm plastid genome annotations, with support for reference-based alignment, quality assessment, and characterisation of intergenic regions. 

## Overview

PAV is designed to process plastid genome assemblies, perform automated annotation using [Chloë](https://github.com/ian-small/Chloe.jl), validate gene annotations against reference sequences, and generate comprehensive reports. The tool includes features for genome linearization, reference-based alignment, detailed quality assessment, and EMBL/ENA template generation.

## Features

- **Automated Annotation**: Uses Chloë for plastid genome annotation
- **Genome Linearization**: Automatically linearizes genomes upstream of a specified gene (defaults to `psbA`)
- **Reference-Based Validation**: Compares annotations against reference genomes from multiple sources (order-specific, default, or custom)
- **Quality Assessment**: Validates gene lengths, identifies internal stop codons, and checks for canonical start and stop codons
- **Alignment Generation**: Creates nucleotide alignments with reference sequences for CDS, rRNA, and tRNA genes
- **Comprehensive Reporting**: Generates detailed reports and statistics
- **EMBL and ENA Template Conversion**: Converts annotated GenBank records to EMBL and produces ENA submission-ready templates
- **Intergenic Region Analysis**: BLAST analysis of intergenic regions for functional characterization

## Installation

### Prerequisites

- Python 3.7+
- [Chloë](https://github.com/ian-small/Chloe.jl) annotation tool
- [MAFFT](https://mafft.cbrc.jp/alignment/software/source.html) (for alignments)
- [TrimAl](https://vicfero.github.io/trimal/index.html) (for backtranslation)
- [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (for querying intergenic regions)

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
  --output_dir /path/to/output
```

### Command Line Options

#### Required Arguments
- `--input_dir`: Directory containing input FASTA files
- `--output_dir`: Directory for output files

#### Optional Arguments
- `--refs_order`: Reference order(s) to use for validation (can specify multiple)
- `--custom_refs_folder`: Custom folder containing reference GenBank files (can be used in addition to --refs_order or default references)
- `--no_alignment`: Skip alignment generation
- `--min_length_percentage`: Minimum gene length percentage (default: 0.8)
- `--max_length_percentage`: Maximum gene length percentage (default: 1.2)
- `--pool`: Number of processes for multiprocessing (default: 1)
- `--threads`: Number of threads per process (default: 1)
- `--linearize_gene`: Gene to use for genome linearization (default: psbA)
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
  --output_dir results/
```

#### Multiple reference orders:
```bash
pav annotate_and_check \
  --input_dir genomes/ \
  --output_dir results/ \
  --refs_order Alismatales Poales Arecales
```

#### Custom reference folder:
```bash
pav annotate_and_check \
  --input_dir genomes/ \
  --output_dir results/ \
  --custom_refs_folder /path/to/custom/references/
```

#### Combined reference sources:
```bash
pav annotate_and_check \
  --input_dir genomes/ \
  --output_dir results/ \
  --refs_order Alismatales \
  --custom_refs_folder /path/to/custom/references/
```

#### Skip alignments:
```bash
pav annotate_and_check \
  --input_dir genomes/ \
  --output_dir results/ \
  --no_alignment
```

#### Custom linearization gene:
```bash
pav annotate_and_check \
  --input_dir genomes/ \
  --output_dir results/ \
  --linearize_gene rbcL
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
- Linearizes genomes upstream of a specified gene (default: psbA)
- Re-annotates linearized sequences

**Linearization Process:**
- The tool identifies the specified gene in the annotated genome
- Linearizes the circular genome by cutting upstream of the gene start position
- Creates a new FASTA file with the linearized sequence
- Re-annotates the linearized sequence to ensure proper feature coordinates

### 2. Reference Validation
- Loads reference sequences from multiple sources (CDS, rRNA, and tRNA)
- Validates gene annotations against reference data
- Checks for gene length outliers, internal stop codons, and non-canonical start and stop codons
- Generates detailed reports

### 3. Alignment Generation
- Creates nucleotide alignments with reference sequences for rRNA and tRNA genes
- Generates protein alignments with backtranslated nucleotide alignments for CDS genes

### 4. EMBL and ENA Template Generation
- Adds locus tags and standardizes features for EMBL
- Builds EMBL template using metadata TSV (see below), ready for submission to ENA

### 6. Intergenic Region Analysis
- Extracts intergenic regions from annotated genomes
- Performs BLAST analysis against reference databases
- Generates comprehensive reports of intergenic region characteristics

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

- **CDS (Coding Sequences)**: Protein-coding genes with translation validation and backtranslation alignment
- **rRNA (Ribosomal RNA)**: Ribosomal RNA genes with direct nucleotide alignment
- **tRNA (Transfer RNA)**: Transfer RNA genes with direct nucleotide alignment

### Reference Sequence Naming
Reference sequences are named using the format: `{Order}_{Family}_{Genus}_{Species}_{GeneName}_{Filename}`
- Multiple copies of the same gene from a single reference file are labeled with `_copy_1`, `_copy_2`, etc.
- This ensures unique identification of each gene copy in alignment files

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

PAV supports multiple sources of reference genomes:

### Built-in References
- **Order-specific references**: Located in `data/order_genomes/` directory, organized by taxonomic order
- **Default references**: Located in `data/reference_genomes_default/` directory

### Custom References
- **Custom reference folder**: Users can provide their own folder of reference GenBank files using `--custom_refs_folder`
- Supports both compressed (.gz) and uncompressed GenBank files
- Can be used alone or in combination with built-in references

### Reference Combination
PAV can combine reference sequences from multiple sources:
- Order-specific references (via `--refs_order`)
- Default references (when no specific references are specified)
- Custom references (via `--custom_refs_folder`)

All reference sources are merged to provide comprehensive validation and alignment data.

### Available Orders
The following taxonomic orders are available in the built-in reference database:
- Acorales, Alismatales, Amborellales, Apiales, Aquifoliales, Arecales, Asparagales, Asterales
- Austrobaileyales, Berberidopsidales, Boraginales, Brassicales, Buxales, Canellales
- Caryophyllales, Celastrales, Ceratophyllales, Chloranthales, Commelinales, Cornales
- Crossosomatales, Cucurbitales, Dilleniales, Dioscoreales, Dipsacales, Ericales
- Fabales, Fagales, Garryales, Gentianales, Geraniales, Huerteales, Icacinales
- Lamiales, Laurales, Liliales, Magnoliales, Malpighiales, Malvales, Metteniusales
- Myrtales, Nymphaeales, Oxalidales, Pandanales, Paracryphiales, Petrosaviales
- Piperales, Poales, Proteales, Ranunculales, Rosales, Santalales, Sapindales
- Saxifragales, Solanales, Trochodendrales, Vitales, Zingiberales, Zygophyllales

See `data/order_genomes/` for the complete list and available reference genomes.

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
6. **Custom reference folder not found**: Verify the path to your custom reference folder exists and contains GenBank files
7. **Linearization gene not found**: If the specified linearization gene is not found in a genome, the original sequence will be used without linearization

## Support

For issues and questions:
- Check the troubleshooting section
- Review the documentation
- Open an issue on GitHub

## Misc links

- ENA EMBL flatfile example: https://ena-docs.readthedocs.io/en/latest/submit/fileprep/flat-file-example.html

## Changelog

### Version 0.0.1
- Initial release
- Genome annotation with Chloë
- Reference-based validation
- Alignment generation
- Quality assessment and reporting
