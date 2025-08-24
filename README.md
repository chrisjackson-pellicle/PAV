            ██████╗      █████╗     ██╗   ██╗
            ██╔══██╗    ██╔══██╗    ██║   ██║
            ██████╔╝    ███████║    ██║   ██║
            ██╔═══╝     ██╔══██║    ╚██╗ ██╔╝
            ██║         ██║  ██║     ╚████╔╝ 
            ╚═╝         ╚═╝  ╚═╝      ╚═══╝  


# Plastid Annotation Validator (PAV)

A tool for annotating and validating angiosperm plastid genome annotations, with support for reference-based alignment, quality assessment, and characterisation of intergenic regions. 

## Overview

PAV is designed to process plastid genome assemblies, perform automated annotation using [Chloë](https://github.com/ian-small/Chloe.jl), validate gene annotations against reference sequences, and generate comprehensive reports. The tool includes features for genome linearisation, reference-based alignment, detailed quality assessment, and EMBL/ENA template generation. PAV supports both single and multi-sequence FASTA files per sample, with the latter allowing for fragmented assemblies.

## Features

- **Automated Annotation**: Uses Chloë for plastid genome annotation
- **Genome Linearisation**: Automatically linearises genomes upstream of a specified gene (defaults to `psbA`) unless genome is specified as 'linear' or the gene is not found
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
- [Chloë](https://github.com/ian-small/Chloe.jl) (annotation tool - not required for `pav check` command)
- [Julia](https://julialang.org/downloads/) (for running Chloë - not required for `pav check` command)
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

   Requires Python packages:

   - biopython>=1.81
   - pandas>=1.5.0
   - tqdm>=4.64.0


   ```bash
   pip install -r requirements.txt
   ```

3. **Install Chloë** (if not already installed):
   ```bash
   # Follow Chloë installation instructions
   # If using a conda environment, ensure chloe and chloe references are in `${CONDA_PREFIX}/bin` or supply the chloe project dir and path to chloe.jl to the `pav annotate_and_check` subcommand
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

#### Full annotation and validation pipeline:
```bash
pav annotate_and_check \
  /path/to/fasta/files \
  /path/to/metadata.tsv \
  --output_directory /path/to/output
```

#### Continue pipeline from annotated GenBank files:
```bash
pav check \
  /path/to/annotated/genbank/files \
  --output_directory /path/to/output
```

See [this caveat](https://github.com/chrisjackson-pellicle/PAV/tree/main?tab=readme-ov-file#using-existing-genbank-files-with-pav-check).

### Command Line Options

#### `pav annotate_and_check` - Full annotation and validation pipeline

```
usage: pav annotate_and_check [-h] [--min_length_percentage FLOAT]
                              [--max_length_percentage FLOAT] [--no_alignment]
                              [--refs_order ORDER [ORDER ...]]
                              [--custom_refs_folder DIR]
                              [--min_intergenic_length INTEGER]
                              [--blast_evalue FLOAT]
                              [--skip_intergenic_analysis]
                              [--debug_intergenic] [--max_blast_hits INTEGER]
                              [--output_directory DIR] [--pool INTEGER]
                              [--threads INTEGER] [--chloe_project_dir DIR]
                              [--chloe_script PATH]
                              [--linearise_gene GENE_NAME] [--run_profiler]
                              DIR TSV

options:
  -h, --help            show this help message and exit

Required input:
  DIR                   Directory containing plastid DNA FASTA files.
  TSV                   TSV file containing sample metadata for EMBL
                        conversion. File should contain columns:
                        input_filename, project_id, locus_tag, genus_species,
                        linear_or_circular. Only input_filename and linear_or_circular 
                        require values; empty optional fields will use defaults.
                        linear_or_circular must be "linear" or "circular".

Optional input:
  --chloe_project_dir DIR, -chloe_proj DIR
                        Path to the chloe project directory to use as
                        --project for Julia. Must be provided together with
                        --chloe_script.
  --chloe_script PATH, -chloe_jl PATH
                        Path to the chloe.jl script. Must be provided together
                        with --chloe_project_dir.
  --linearise_gene GENE_NAME, -linearise_gene GENE_NAME
                        Gene to use for genome linearisation. Default is: psbA
```

#### `pav check` - Continue pipeline from annotated GenBank files

```
usage: pav check [-h] [--min_length_percentage FLOAT]
                 [--max_length_percentage FLOAT] [--no_alignment]
                 [--refs_order ORDER [ORDER ...]]
                 [--custom_refs_folder DIR]
                 [--min_intergenic_length INTEGER]
                 [--blast_evalue FLOAT]
                 [--skip_intergenic_analysis]
                 [--debug_intergenic] [--max_blast_hits INTEGER]
                 [--output_directory DIR] [--pool INTEGER]
                 [--threads INTEGER] [--metadata_tsv TSV]
                 [--run_profiler]
                 DIR

options:
  -h, --help            show this help message and exit

Required input:
  DIR                   Directory containing annotated GenBank files (.gb, .gbk, .gb.gz, .gbk.gz).

Optional input:
  --metadata_tsv TSV    TSV file containing sample metadata for EMBL conversion.
                        If not provided, samples will be processed without metadata.
```

#### Options common to both `pav annotate_and_check` and `pav check`

```
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
  --refs_order ORDER [ORDER ...], -refs_ord ORDER [ORDER ...]
                        Order(s) to use for reference genes. Can be specified
                        multiple times. Default is: []
  --custom_refs_folder DIR, -custom_refs DIR
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
  --output_directory DIR, -out_dir DIR
                        Output directory for the subcommand. Default is:
                        output_directory
  --pool INTEGER, -p INTEGER
                        The number of CPUs to use for the subcommand. Default
                        is: 1
  --threads INTEGER, -t INTEGER
                        The number of threads to use for the subcommand.
                        Default is: 1
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

#### Continue from annotated GenBank files:
```bash
pav check \
  annotated_genomes/
```

#### Continue from annotated GenBank files with metadata:
```bash
pav check \
  annotated_genomes/ \
  --metadata_tsv metadata.tsv
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

#### Custom linearization gene(s):
```bash
# Single gene (default behavior)
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --linearise_gene psbA

# Multiple genes (tried in order until one is found)
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --linearise_gene psbA rbcL matK

# Mixed case (handles gene name mapping)
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --linearise_gene ATPA psbA rbcL
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

### For `pav annotate_and_check`:
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
│       ├── <prefix>.round1.chloe.gbk            # Original annotation (preserved)
│       ├── <prefix>.round1.chloe.gff            # Original GFF (preserved)
│       ├── <prefix>.round2.chloe.gbk            # Re-annotated after linearization
│       ├── <prefix>.round2.chloe.gff            # Re-annotated after linearization
│       └── <prefix>.round2.fasta                # Linearised sequence
│       └── <prefix>.round2.chloe_intergenic_debug.fasta   # Optional (when --debug_intergenic)
│       └── <prefix>_seq001_<seqname>.fasta      # Individual sequences (multi-sequence files only)
│       └── <prefix>_seq001_<seqname>.round1.chloe.gbk  # Individual sequence annotations
│       └── <prefix>_seq001_<seqname>.round1.chloe.gff  # Individual sequence GFF files
├── 02_embl_files/
│   ├── <sample_name>.embl                        # EMBL format
│   └── ena_submission_embl_files/
│       └── <sample_name>.ena.embl                # ENA template (derived from EMBL)
├── 03_alignments_with_refs/
│   ├── 01_per_sample_alignments/
│   │   └── <sample_name>/
│   │       ├── <sample>_<gene>_alignment.fasta          # Per-sample CDS alignment
│   │       ├── <sample>_<gene>_rRNA_alignment.fasta     # Per-sample rRNA alignment
│   │       └── <sample>_<gene>_tRNA_alignment.fasta     # Per-sample tRNA alignment
│   └── 02_per_gene_alignments/
│       └── <gene>_all_samples_alignment.fasta           # All samples combined (nucleotide)
├── 04_intergenic_analysis/
    ├── <sample_name>_intergenic_blast_results.tsv       # Per-sample BLAST results
    └── combined_intergenic_blast_results.tsv            # Combined BLAST results
```

### For `pav check`:
```
output_dir/
├── 00_logs_and_reports/
│   ├── logs/
│   │   └── check_<timestamp>.log                         # Runtime log
│   └── reports/
│       ├── all_samples_gene_validation_report.tsv       # Combined gene length/translation checks
│       └── <sample_name>_gene_validation_report.tsv     # Per-sample reports
├── 01_gbk_split/                                        # Split multi-record GenBank files
│   └── <sample_name>_seq001_<seqname>.gb                # Individual sequence files from multi-record GenBank files
├── 02_embl_files/
│   ├── <sample_name>.embl                               # EMBL format
│   └── ena_submission_embl_files/
│       └── <sample_name>.ena.embl                       # ENA template (derived from EMBL)
├── 03_alignments_with_refs/
│   ├── 01_per_sample_alignments/
│   │   └── <sample_name>/
│   │       ├── <sample>_<gene>_alignment.fasta          # Per-sample CDS alignment
│   │       ├── <sample>_<gene>_rRNA_alignment.fasta     # Per-sample rRNA alignment
│   │       └── <sample>_<gene>_tRNA_alignment.fasta     # Per-sample tRNA alignment
│   └── 02_per_gene_alignments/
│       └── <gene>_all_samples_alignment.fasta           # All samples combined (nucleotide)
├── 04_intergenic_analysis/
    ├── <sample_name>_intergenic_blast_results.tsv       # Per-sample BLAST results
    └── combined_intergenic_blast_results.tsv            # Combined BLAST results
```

## Workflow

### `pav annotate_and_check` - Full Pipeline

#### 1. Genome Annotation
- Processes input FASTA files using Chloë
- Supports both single and multi-sequence FASTA files per sample
- For multi-sequence FASTA files, each sequence is processed separately with Chloë
- Performs initial annotation on original sequences
- Linearises genomes upstream of specified gene(s) (default: `psbA`), unless sample is recorded as `linear` in metadata. Multiple genes can be specified and will be tried in order until one is found. If no specified genes are found in the sequence, no linearisation occurs
- Re-annotates linearised sequences

#### 2. Reference Validation
- Loads reference sequences from multiple sources (CDS, rRNA, and tRNA)
- Validates gene annotations against reference data
- Checks for gene length outliers, internal stop codons, and non-canonical start and stop codons
- Generates detailed reports

#### 3. Alignment Generation
- Creates nucleotide alignments with reference sequences for rRNA and tRNA genes
- Generates codon-aware nucleotide alignments for CDS genes

#### 4. EMBL and ENA Template Generation
- Adds locus tags and standardizes features for EMBL
- Builds EMBL templates using metadata TSV (see below), ready for submission to ENA

#### 5. Intergenic Region Analysis
- Extracts intergenic regions from annotated genomes
- Performs BLAST analysis against a reference database of coding regions
- Identifies potential functional intergenic regions
- Generates comprehensive reports of intergenic region characteristics

### `pav check` - Continue from Annotated GenBank Files

#### 1. Load Annotated GenBank Files
- Processes existing annotated GenBank files (.gb, .gbk, .gb.gz, .gbk.gz)
- Automatically splits multi-record GenBank files into individual sequence files
- Supports both single and multi-record files per sample

#### 2. Reference Validation
- Same validation process as the full pipeline
- Validates gene annotations against reference data
- Checks for gene length outliers, internal stop codons, and non-canonical start and stop codons

#### 3. Alignment Generation
- Same alignment process as the full pipeline
- Creates nucleotide alignments with reference sequences

#### 4. EMBL and ENA Template Generation
- Same EMBL conversion process as the full pipeline
- Requires metadata TSV file if provided via `--metadata_tsv`

#### 5. Intergenic Region Analysis
- Same intergenic analysis as the full pipeline
- Extracts and analyzes intergenic regions

## Linearization

PAV automatically linearises plastid genomes upstream of specified gene(s) to ensure consistent annotation and to avoid issues with genes that span the circular genome boundary.

## Gene naming conventions

PAV follows the gene naming conventions used by Chloë, e.g.

- `pafI` rather than `ycf3`
- `pbf1` rather than `psbN`
- `ndhK` rather than `psbG`
- etc.

This means that if you run `pav check` on a genome with the annotation `ycf3`, the corresponding alignment produced by PAV will be named e.g. `pafI_all_samples_alignment.fasta`.

See the full list of synonyms [here](https://github.com/chrisjackson-pellicle/PAV/blob/main/plastid_annotation_validator/data/gene_synonyms.txt).

## Metadata TSV Format

Provide a tab-separated file with the following columns (include this header in the file):

- `input_filename` (required)
- `project_id` (optional)
- `locus_tag` (optional)
- `genus_species` (optional)
- `linear_or_circular` (required, must be either `linear` or `circular`)

Example [`metadata.tsv`](https://github.com/chrisjackson-pellicle/PAV/blob/main/plastid_annotation_validator/data/metadata.tsv):

| input_filename | project_id | locus_tag | genus_species | linear_or_circular |
|----------------|------------|-----------|---------------|--------------------|
| sample1.fasta  | PRJEB12345 | ABC       | Arabidopsis thaliana | circular |
| sample2.fasta  | PRJEB98765 | XYZ       | Oryza sativa | linear |
| sample3.fasta  |            |           |               | circular |
| sample4.fasta  | PRJEB11111 |           | Zea mays | circular |


Notes:
- **Required fields**: `input_filename` and `linear_or_circular` must have values when running `pav annotate_and_check`
- **Optional fields**: `project_id`, `locus_tag`, and `genus_species` can be empty
- **Default values** for empty optional fields:
  - `project_id` → `'UNKNOWN_PROJECT'`
  - `locus_tag` → `'DEFAULT_TAG'`
  - `genus_species` → `'Unknown species'`
  - `linear_or_circular` → `'unknown'` (only allowed when running `pav check`)
- **File format**: The metadata file must be a tab-separated file with the above columns
- **File path**: The path to the metadata file is a required positional argument for `pav annotate_and_check`, but an optional parameter passed via `--metadata_tsv` for `pav check`
- The topology is used in the ENA `ID` line and must match one of: `linear`, `circular`
- **For `pav annotate_and_check`**: All samples must be present in the metadata file
- **For `pav check`**: Metadata is optional - if not provided, samples will be processed without metadata
- Samples marked as `linear` will **not be relinearised** at the `--linearise_gene` position, preserving their original linear structure

## Intergenic regions reports

Note that the database of coding regions used in the intergenic analyses is derived from the genomes in the [`/data/order_genomes`](https://github.com/chrisjackson-pellicle/PAV/tree/main/plastid_annotation_validator/data/order_genomes) folder. These genomes can have mis-annotations, leading to **spurious detection of genes in the intergenic regions of your genomes**. For example, if an `rrn23` annotation in the reference database has a 5' termini that was incorrectly extended in to upstream intergenic sequence (e.g. `NC_036304`), `rrn23` may be reported in the homologous intergenic region of your genome. This is a limitation of the reference database, and will be improved in future releases. 

## Using existing Genbank files with `pav check`

PAV can be used to check the annotations of existing plastid genome GenBank files. To use this feature, provide a folder of GenBank files as input to the `pav check` command. 

However, some annotation programs (including Chloe) attempt to annotate genes that might cross/bridge the ends of a linear fasta input sequence (e.g. the 5' end of `rps19` will be annotated at the 3' end of the input fasta sequence, and the 3' end of `rps19` will be annotated at the 5' end of the fasta sequence). In the case of Chloë (as of 23rd August 2025), if an incomplete fasta sequence (e.g. a sequence that can not be circularised, and has a partial `rps19` annotation at the 3' end of the fasta sequence, whereas the 5' end begins with intergenic sequence) this can result in a partially incorrect `rps19` annotation. This occurs because Chloe identifies the CDS for a putative gene by finding a corresponding ORF; if the 5' end of the fasta sequence (intergenic sequence in this example) allows the ORF of the correct `rps19` region at the 3' end of the fasta sequence to be extended before reaching an in-frame stop codon, the resulting ORF (containing intergenic sequence) will be returned as the `rps19` CDS.

Even if the `rps19` annotation is correct (i.e. the input fasta sequence can be circularised, but happens to be linearised halfway through the `rps19` CDS), the outcome is an annotation interval that extends beyond the end length coordinate of the input sequence (e.g. `rps19` is annotated with an interval of `156176..156472`, whereas the input sequence is only 156,374 bp). This is an issue with the current implementation of PAV, which uses the Biopython SeqIO Genbank parser to extract CDS sequences (i.e. using `feature.location.extract(record.seq)`). The Biopython Genbank parser does not handle the interval `156176..156472` correctly, and instead returns a sequence corresonding to the interval `156176..156374` (i.e. the sequence is truncated to the 3' end of the fasta sequence). In these scenarios, PAV produces gene length warnings in the report files, such as:

```
['Gene rps19 length mismatch: Expected 297 != Extracted 199 for copy 1 - check this gene manually!', 'CDS length 199 is not a multiple of 3 for copy 1', "Stop codon is TGG, expected one of ['TAA', 'TAG', 'TGA'] for copy 1"]
```

When annotating fasta sequences using `pav annotate_and_check`, PAV tries to avoid this issue using a two-pass approach:
 1) The input sequence is annotated once with Chloe
 2) The location of a specified gene is recovered (`psbA` by default, can be changed with `--linearise_gene <gene_name>`) and the input fasta is linearised just upstream.
 3) The input sequence is annotated again with Chloe, this time using the linearised sequence as input. Consequently, no annotation should cross/bridge the ends of the annotated fasta sequence. 

Note that if the specified gene for linearisation can't be found, the fasta sequence is processed as is. So, if a fragmented assembly with multiple fasta records is used as input, `psbA` may be found (and used for linearisation) in one sequence, but the others will be processed as is (i.e. no linearisation is applied), potentially leading to the issue described above. 


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
- **FASTA**: Genome assemblies (.fasta, .fa, .fas) - for `pav annotate_and_check`
- **GenBank**: Annotated genomes (.gb, .gbk, .gb.gz, .gbk.gz) - for `pav check`

### Output
- **GenBank**: Annotated genomes in GenBank format (from Chloë)
- **GFF**: Gene feature format files (from Chloë)
- **EMBL**: European Molecular Biology Laboratory format
- **FASTA**: Aligned sequences and linearised genomes

## Reference Data

PAV supports multiple sources of reference genomes:

### Built-in References
- **Order-specific references**: Located in [`data/order_genomes`](https://github.com/chrisjackson-pellicle/PAV/tree/main/plastid_annotation_validator/data/order_genomes) directory, organized by taxonomic order
- **Default references**: Located in [`data/reference_genomes_default`](https://github.com/chrisjackson-pellicle/PAV/tree/main/plastid_annotation_validator/data/reference_genomes_default) directory

### Custom References
- **Custom reference folder**: Users can provide their own folder of reference GenBank files using `--custom_refs_folder`
- Supports both compressed (`.gz`) and uncompressed GenBank files
- Can be used alone or in combination with built-in references

### Reference Combination
PAV can combine reference sequences from multiple sources:
- Order-specific references (via `--refs_order`)
- Custom references (via `--custom_refs_folder`)
- Default references (when no specific references are specified)

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

See [`data/order_genomes`](https://github.com/chrisjackson-pellicle/PAV/tree/main/plastid_annotation_validator/data/order_genomes) for the complete list and available reference genomes.

## Configuration

### Median Lengths
Reference gene lengths are stored in [`data/plDNA_genes_median_lengths.csv`](https://github.com/chrisjackson-pellicle/PAV/blob/main/plastid_annotation_validator/data/plDNA_genes_median_lengths.csv) for validation.

## Troubleshooting

### Common Issues

1. **Chloë not found**: Ensure Chloë is properly installed and accessible
2. **MAFFT/TrimAl errors**: Check that external tools are installed and in your $PATH
3. **BLAST+ not found**: Install BLAST+ (`blastn`, `makeblastdb`) and ensure they are in your $PATH
4. **Metadata TSV columns**: Ensure the TSV has the required columns listed above
5. **Memory issues**: Reduce `--pool` and `--threads` parameters
6. **Custom reference folder not found**: Verify the path to your custom reference folder exists and contains GenBank files
7. **Linearization gene(s) not found**: If none of the specified linearization genes are found in a genome, the original sequence will be used without linearization. Genes are validated against the gene_synonyms.txt file to ensure they are legitimate plastid genes
8. **Gzipped GenBank files**: PAV supports both compressed (.gz) and uncompressed GenBank files as input for the `check` subcommand

## Support

For issues and questions:
- Check the Troubleshooting section
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
