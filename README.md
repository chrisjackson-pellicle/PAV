            ██████╗      █████╗     ██╗   ██╗
            ██╔══██╗    ██╔══██╗    ██║   ██║
            ██████╔╝    ███████║    ██║   ██║
            ██╔═══╝     ██╔══██║    ╚██╗ ██╔╝
            ██║         ██║  ██║     ╚████╔╝ 
            ╚═╝         ╚═╝  ╚═╝      ╚═══╝  


# Plastid Annotation Validator (PAV)

A tool for annotating and validating angiosperm plastid genome annotations, with support for reference-based alignment, quality assessment, and characterisation of intergenic regions. 

## Overview

PAV is designed to process plastid genome assemblies, perform automated annotation using [Chloë](https://github.com/ian-small/Chloe.jl), validate gene annotations against reference sequences, and generate comprehensive reports. The tool includes features for genome linearization, reference-based alignment, detailed quality assessment, and EMBL/ENA template generation. PAV supports both single and multi-sequence FASTA files, automatically detecting and processing each sequence individually.

## Features

- **Automated Annotation**: Uses Chloë for plastid genome annotation
- **Genome Linearisation**: Automatically linearises genomes upstream of a specified gene (defaults to `psbA`)
- **Reference-Based Validation**: Compares annotations against reference genomes from multiple sources (order-specific, default, or custom)
- **Quality Assessment**: Validates gene lengths, identifies internal stop codons, and checks for canonical start and stop codons
- **Alignment Generation**: Creates nucleotide alignments with reference sequences for CDS, rRNA, and tRNA genes
- **Comprehensive Reporting**: Generates detailed reports and statistics
- **EMBL and ENA Template Conversion**: Converts annotated GenBank records to EMBL and produces ENA submission-ready templates
- **Intergenic Region Analysis**: BLAST analysis of intergenic regions for functional characterization
- **Multi-Sequence Support**: Automatic detection and processing of multi-sequence FASTA files

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
   pip install .
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
  /path/to/fasta/files \
  /path/to/metadata.tsv \
  --output_directory /path/to/output
```

### Command Line Options

```
usage: pav annotate_and_check [-h] [--min_length_percentage FLOAT]
                              [--max_length_percentage FLOAT] [--no_alignment]
                              [--refs_order REFS_ORDER [REFS_ORDER ...]]
                              [--custom_refs_folder CUSTOM_REFS_FOLDER]
                              [--min_intergenic_length INTEGER]
                              [--blast_evalue FLOAT]
                              [--skip_intergenic_analysis]
                              [--debug_intergenic] [--max_blast_hits INTEGER]
                              [--output_directory output_directory]
                              [--pool INTEGER] [--threads INTEGER]
                              [--chloe_project_dir PATH] [--chloe_script PATH]
                              [--linearize_gene GENE_NAME] [--run_profiler]
                              DIR TSV

options:
  -h, --help            show this help message and exit

Required input:
  DIR                   Directory containing plastid DNA FASTA files (supports both single and multi-sequence files).
  TSV                   TSV file containing sample metadata for EMBL
                        conversion. Required file should contain columns:
                        fasta_filename, project_id, locus_tag, genus_species,
                        linear_or_circular. ALL samples must be listed. Only
                        fasta_filename and linear_or_circular require values;
                        empty optional fields will use defaults.
                        linear_or_circular must be "linear" or "circular".

Gene length warnings:
  --min_length_percentage FLOAT, -min_len FLOAT
                        Minimum length percentage of the gene median length
                        for a warning to be issued. Default is: 0.8
  --max_length_percentage FLOAT, -max_len FLOAT
                        Maximum length percentage of the gene median length
                        for a warning to be issued. Default is: 1.2

Alignment with reference genes:
  --no_alignment, -no_align
                        Do not align annotated genes with reference genes.
                        Default is: False
  --refs_order REFS_ORDER [REFS_ORDER ...], -refs_ord REFS_ORDER [REFS_ORDER ...]
                        Order(s) to use for reference genes. Can be specified
                        multiple times. Default is: []
  --custom_refs_folder CUSTOM_REFS_FOLDER, -custom_refs CUSTOM_REFS_FOLDER
                        Custom folder containing reference GenBank files. Can
                        be used in addition to --refs_order or default
                        references.

Intergenic region analysis:
  --min_intergenic_length INTEGER, -min_ig_len INTEGER
                        Minimum length of intergenic region to analyze.
                        Default is: 0
  --blast_evalue FLOAT, -evalue FLOAT
                        BLAST E-value threshold for intergenic region
                        analysis. Default is: 1e-10
  --skip_intergenic_analysis, -skip_ig
                        Skip intergenic region analysis. Default is: False
  --debug_intergenic, -debug_ig
                        Write intergenic regions to FASTA files for debugging.
                        Default is: False
  --max_blast_hits INTEGER, -max_hits INTEGER
                        Maximum number of BLAST hits to retain per intergenic
                        region. Default is: 1

General pipeline options:
  --output_directory output_directory, -out_dir output_directory
                        Output directory for the subcommand. Default is:
                        output_directory
  --pool INTEGER        The number of CPUs to use for the subcommand. Default
                        is: 1
  --threads INTEGER     The number of threads to use for the subcommand.
                        Default is: 1
  --chloe_project_dir PATH, -chloe_proj PATH
                        Path to the chloe project directory to use as
                        --project for Julia. Must be provided together with
                        --chloe_script.
  --chloe_script PATH, -chloe_jl PATH
                        Path to the chloe.jl script. Must be provided together
                        with --chloe_project_dir.
  --linearize_gene GENE_NAME, -linearize_gene GENE_NAME
                        Gene to use for genome linearization. Default is: psbA
  --run_profiler        If supplied, run the subcommand using cProfile. Saves
                        a *.csv file of results. Default is: False
```

### Example Commands

#### Basic annotation and validation:
```bash
pav annotate_and_check \
  genomes/ \
  metadata.tsv
```

#### Multiple reference orders:
```bash
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --refs_order Alismatales Poales Arecales
```

#### Custom reference folder:
```bash
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --custom_refs_folder /path/to/custom/references/
```

#### Combined reference sources:
```bash
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --refs_order Alismatales \
  --custom_refs_folder /path/to/custom/references/
```

#### Skip alignments:
```bash
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --no_alignment
```

#### Custom linearization gene:
```bash
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --linearize_gene rbcL
```

#### High-performance processing:
```bash
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --output_directory results/ \
  --refs_order Alismatales \
  --min_length_percentage 0.9 \
  --max_length_percentage 1.1 \
  --pool 8 \
  --threads 4
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
│       └── <prefix>_seq001_<seqname>.fasta      # Individual sequences (multi-sequence files only)
│       └── <prefix>_seq001_<seqname>.chloe.gbk  # Individual sequence annotations
│       └── <prefix>_seq001_<seqname>.chloe.gff  # Individual sequence GFF files
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
- **Supports both single and multi-sequence FASTA files**
- For multi-sequence FASTA files, each sequence is processed separately with Chloë
- Performs initial annotation on original sequences
- Linearizes genomes upstream of a specified gene (default: psbA), unless sample is recorded as `linear` in metadata
- Re-annotates linearized sequences

### 2. Reference Validation
- Loads reference sequences from multiple sources (CDS, rRNA, and tRNA)
- Validates gene annotations against reference data
- Checks for gene length outliers, internal stop codons, and non-canonical start and stop codons
- Generates detailed reports

### 3. Alignment Generation
- Creates nucleotide alignments with reference sequences for rRNA and tRNA genes
- Generates codon-aware nucleotide alignments for CDS genes

### 4. EMBL and ENA Template Generation
- Adds locus tags and standardizes features for EMBL
- Builds EMBL templates using metadata TSV (see below), ready for submission to ENA

### 5. Intergenic Region Analysis
- Extracts intergenic regions from annotated genomes
- Performs BLAST analysis against reference databases
- Generates comprehensive reports of intergenic region characteristics

## Metadata TSV Format

Provide a tab-separated file with the following required columns (include this header in the file):

- `fasta_filename`
- `project_id`
- `locus_tag`
- `genus_species`
- `linear_or_circular` (must be either `linear` or `circular`)

Example (`metadata.tsv`):

| fasta_filename | project_id | locus_tag | genus_species | linear_or_circular |
|----------------|------------|-----------|---------------|--------------------|
| sample1.fasta  | PRJEB12345 | ABC       | Arabidopsis thaliana | circular |
| sample2.fasta  | PRJEB98765 | XYZ       | Oryza sativa | linear |
| sample3.fasta  |            |           |               | circular |


Notes:
- All columns above must be present in the header, but only `fasta_filename` and `linear_or_circular` require values
- Empty optional fields will use these default values:
  - `project_id` → `'UNKNOWN_PROJECT'`
  - `locus_tag` → `'DEFAULT_TAG'`
  - `genus_species` → `'Unknown species'`
- The topology is used in the ENA `ID` line and must match one of: `linear`, `circular`
- **ALL samples must be present in the metadata file** - PAV will fail if any samples are missing
- Samples marked as `linear` will **not be re-linearized** at the `--linearize_gene` position, preserving their original linear structure

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

## Misc. links

- ENA EMBL flatfile example: https://ena-docs.readthedocs.io/en/latest/submit/fileprep/flat-file-example.html

## Changelog

### Version 0.0.1
- Initial release
- Genome annotation with Chloë
- Reference-based validation
- Alignment generation
- Quality assessment and reporting
