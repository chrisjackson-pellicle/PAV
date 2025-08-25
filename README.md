            ██████╗      █████╗     ██╗   ██╗
            ██╔══██╗    ██╔══██╗    ██║   ██║
            ██████╔╝    ███████║    ██║   ██║
            ██╔═══╝     ██╔══██║    ╚██╗ ██╔╝
            ██║         ██║  ██║     ╚████╔╝ 
            ╚═╝         ╚═╝  ╚═╝      ╚═══╝  


# Plastid Annotation Validator (PAV)

A tool for annotating and validating angiosperm plastid genome annotations, with support for reference-based alignment, quality assessment, and characterisation of intergenic regions. 

## Overview

PAV processes plastid genome assemblies, performs automated annotation using [Chloë](https://github.com/ian-small/Chloe.jl), validates gene annotations against reference sequences, and generates comprehensive reports. The tool includes features for genome linearisation, reference-based alignment, detailed quality assessment, and EMBL/ENA template generation. PAV supports both single and multi-sequence FASTA files per sample, with the latter allowing for fragmented assemblies.

## Features

- **Automated annotation**: uses Chloë for plastid genome annotation
- **Genome linearisation**: automatically linearises genomes upstream of a specified gene (defaults to `psbA`) unless genome is specified as 'linear' or the gene is not found
- **Quality assessment**: validates gene lengths, identifies internal stop codons, and checks for canonical start and stop codons
- **Reference-based validation**: compares annotations against reference genomes from multiple sources (order-specific, default, or custom)
- **Alignment generation**: creates nucleotide alignments with reference sequences for CDS, rRNA, and tRNA genes
- **Comprehensive reporting**: generates detailed reports and statistics
- **EMBL and ENA template conversion**: converts annotated GenBank records to EMBL and produces ENA submission-ready templates
- **Intergenic region analysis**: BLAST analysis of intergenic regions for functional characterization, with custom database option
- **Multi-sequence support**: automatic detection and processing of multi-sequence FASTA files

## Installation

### Prerequisites

- Python 3.7+
- [Chloë](https://github.com/ian-small/Chloe.jl) (annotation tool - not required for `pav check` command)
- [Julia](https://julialang.org/downloads/) (for running Chloë - not required for `pav check` command)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/source.html) (for alignments)
- [TrimAl](https://vicfero.github.io/trimal/index.html) (for backtranslation)
- [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (for querying intergenic regions)

### Setup

1. **Clone the PAV repository**:
   ```bash
   git clone https://github.com/chrisjackson-pellicle/PAV.git
   cd PAV
   pip install .
   ```

2. **Install Python dependencies**:
  <br/><br/>
   Requires Python packages:

   - biopython>=1.81
   - pandas>=1.5.0
   - tqdm>=4.64.0

   ```bash
   pip install -r requirements.txt
   ```

3. **Install Chloë**:
  <br/><br/>
   Follow Chloë installation instructions [here](https://github.com/ian-small/Chloe.jl?tab=readme-ov-file#installation)
  <br/><br/>
   If using a conda environment, ensure Chloe and Chloe references are in `${CONDA_PREFIX}/bin` directory or supply the Chloe project dir and path to `Chloe.jl` to the `pav annotate_and_check` subcommand using the options `--chloe_project_dir` and `--chloe_script`, respectively. 


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

### Basic usage

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

#### Create custom BLAST database for intergenic analysis:
```bash
pav db_for_intgen \
  /path/to/genbank/files \
  --output_directory /path/to/output
```

### Command line options

#### `pav annotate_and_check` - full annotation and validation pipeline

```
usage: pav annotate_and_check [-h] [--min_length_percentage FLOAT]
                              [--max_length_percentage FLOAT] [--no_alignment]
                              [--refs_order ORDER [ORDER ...]]
                              [--custom_refs_folder DIR]
                              [--min_intergenic_length INTEGER]
                              [--blast_evalue FLOAT]
                              [--skip_intergenic_analysis]
                              [--debug_intergenic] [--max_blast_hits INTEGER]
                              [--custom_blast_db PATH]
                              [--output_directory DIR] [--pool INTEGER]
                              [--threads INTEGER] [--chloe_project_dir DIR]
                              [--chloe_script PATH]
                              [--linearise_gene GENE_NAME [GENE_NAME ...]] 
                              [--run_profiler]
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
  --linearise_gene GENE_NAME [GENE_NAME ...]
                        Gene(s) to use for genome linearisation. Can specify multiple genes.
                        Genes will be tried in order until one is found. Default is: ['psbA']
```

#### `pav check` - continue pipeline from annotated GenBank files

```
usage: pav check [-h] [--min_length_percentage FLOAT]
                 [--max_length_percentage FLOAT] [--no_alignment]
                 [--refs_order ORDER [ORDER ...]]
                 [--custom_refs_folder DIR]
                 [--min_intergenic_length INTEGER]
                 [--blast_evalue FLOAT]
                 [--skip_intergenic_analysis]
                 [--debug_intergenic] [--max_blast_hits INTEGER]
                 [--custom_blast_db PATH]
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
  --custom_blast_db PATH, -blast_db PATH
                        Custom BLAST database path or directory for intergenic
                        region analysis. Can be a direct path to the database
                        or a directory containing BLAST database files. If not
                        provided, uses the default order_genomes_blastdb.

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

### Example commands

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

#### Custom BLAST database for intergenic analysis:
```bash
# Using a directory containing BLAST database files
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --custom_blast_db /path/to/blastdb/directory/

# Using a direct path to the BLAST database
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --custom_blast_db /path/to/blastdb/directory/custom_db
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

#### Create custom BLAST database for intergenic analysis:
```bash
# Basic usage
pav db_for_intgen \
  genbank_files/ \
  --output_directory blastdb_output/

# With custom minimum length filter
pav db_for_intgen \
  genbank_files/ \
  --output_directory blastdb_output/ \
  --min_length 50
```

#### Using custom BLAST database with other subcommands:
```bash
# First, create a custom BLAST database
pav db_for_intgen \
  reference_genbank_files/ \
  --output_directory custom_blastdb/

# Then use it for intergenic analysis
pav annotate_and_check \
  genomes/ \
  metadata.tsv \
  --custom_blast_db custom_blastdb/01_intgen_db/

# Or with the check subcommand
pav check \
  annotated_genomes/ \
  --custom_blast_db custom_blastdb/01_intgen_db/
```

## Output structure

PAV generates a structured output directory:

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

### For `pav db_for_intgen`:
```
output_dir/
├── 00_logs_and_reports/
│   ├── logs/
│   │   └── db_for_intgen_<timestamp>.log                # Runtime log
└── 01_intgen_db/
    ├── intgen_db.nhr                                    # BLAST database header file
    ├── intgen_db.nin                                    # BLAST database index file
    ├── intgen_db.nsq                                    # BLAST database sequence file
    ├── intgen_db.ndb                                    # BLAST database index file
    ├── intgen_db.not                                    # BLAST database note file
    ├── intgen_db.ntf                                    # BLAST database taxonomy file
    └── intgen_db.nto                                    # BLAST database taxonomy offset file
```

## Workflow

### `pav annotate_and_check` - full pipeline

#### 1. Genome annotation
- Processes input FASTA files using Chloë
- Supports both single and multi-sequence FASTA files per sample
- For multi-sequence FASTA files, each sequence is processed separately with Chloë
- Performs initial annotation on original sequences
- Linearises genomes upstream of specified gene(s) (default: `psbA`), unless sample is recorded as `linear` in metadata. Multiple genes can be specified and will be tried in order until one is found. If no specified genes are found in the sequence, no linearisation occurs
- Re-annotates linearised sequences

#### 2. Annotation validation
- Validates gene annotations against reference data
- Checks for gene length outliers, internal stop codons, and non-canonical start and stop codons
- Generates detailed reports

#### 3. Alignment generation
- Creates nucleotide alignments with reference sequences for rRNA and tRNA genes
- Generates codon-aware nucleotide alignments for CDS genes

#### 4. EMBL and ENA template generation
- Adds locus tags and standardizes features for EMBL
- Builds EMBL templates using metadata TSV (see below), ready for submission to ENA

#### 5. Intergenic region analysis
- Extracts intergenic regions from annotated genomes
- Performs BLAST analysis against a reference database of coding regions, or a custom database
- Identifies potential functional intergenic regions
- Generates comprehensive reports of intergenic region characteristics

### `pav check` - continue from annotated GenBank Files

#### 1. Load annotated GenBank files
- Processes existing GenBank files with annotations (.gb, .gbk, .gb.gz, .gbk.gz)
- Supports both single and multi-record files per sample
- Automatically splits multi-record GenBank files into individual sequence files

#### 2. Reference validation
- Same validation process as the full pipeline
- Validates gene annotations against reference data
- Checks for gene length outliers, internal stop codons, and non-canonical start and stop codons

#### 3. Alignment generation
- Same alignment process as the full pipeline
- Creates nucleotide alignments with reference sequences

#### 4. EMBL and ENA template generation
- Same EMBL conversion process as the full pipeline
- Metadata TSV file can optionally be provided via `--metadata_tsv`

#### 5. Intergenic region analysis
- Same intergenic analysis as the full pipeline
- Extracts and analyzes intergenic regions

### `pav db_for_intgen` - create custom BLAST database

#### 1. GenBank file processing
- Processes GenBank files (.gb, .gbk, .gb.gz, .gbk.gz) from the input directory
- Extracts gene sequences and features from each GenBank file
- Excludes regions: inverted repeats, LSC, SSC, repeat_region, misc_feature, intron, gene, and source features

#### 2. Feature extraction and filtering
- Extracts CDS, rRNA, tRNA, and other functional features
- Filters features by minimum length (default: 10 bp, configurable via `--min_length`)

#### 3. Sequence ID generation
- Creates standardized sequence IDs in the format: `Order_Family_Genus_Species_accession_feature_type_[gene_name]`
- Extracts taxonomic information directly from GenBank file annotations
- Handles cases where taxonomic information is incomplete

#### 4. BLAST database creation
- Writes extracted features to a temporary FASTA file
- Creates a nucleotide BLAST database using `makeblastdb`

## Linearization

PAV automatically linearises plastid genomes upstream of specified gene(s) to ensure consistent annotation and to avoid issues with genes that span the circular genome boundary (see below).

## Gene naming conventions

PAV follows the gene naming conventions used by Chloë, e.g.

- `pafI` rather than `ycf3`
- `pbf1` rather than `psbN`
- `ndhK` rather than `psbG`
- etc.

This means that if you run `pav check` on a genome with the annotation `ycf3`, the corresponding alignment produced by PAV will be named e.g. `pafI_all_samples_alignment.fasta`.

See the full list of synonyms used for mapping gene names [here](https://github.com/chrisjackson-pellicle/PAV/blob/main/plastid_annotation_validator/data/gene_synonyms.txt).

## Metadata TSV format

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
- **File format**: the metadata file must be a tab-separated file with the above columns
- **File path**: the path to the metadata file is a required positional argument for `pav annotate_and_check`, but an optional parameter passed via `--metadata_tsv` for `pav check`
- The topology is used in the ENA `ID` line and must match one of: `linear`, `circular`
- **For `pav annotate_and_check`**: all samples must be present in the metadata file
- **For `pav check`**: metadata is optional - if not provided, samples will be processed without metadata
- Samples marked as `linear` will **not be relinearised** at the `--linearise_gene` position, preserving their original linear structure

## Intergenic regions reports

Note that the database of coding regions used in the intergenic analyses is derived from the genomes in the [`/data/order_genomes`](https://github.com/chrisjackson-pellicle/PAV/tree/main/plastid_annotation_validator/data/order_genomes) folder. These genomes can have mis-annotations, leading to **spurious detection of genes in the intergenic regions of your genomes**. For example, if an `rrn23` annotation in the reference database has a mis-annotated5' termini that was incorrectly extended in to upstream intergenic sequence (e.g. `NC_036304`), `rrn23` may be reported in the homologous intergenic region of your genome. This is a limitation of the reference database, and will be improved in future releases. 

In the meantime, users can create their own database of coding regions by running `pav db_for_intgen` on a set of trusted plastid genome GenBank files, and then using the resulting database for intergenic analyses in `pav annotate_and_check` or `pav check` via the `--custom_blast_db` option.   

## Using existing GenBank files with `pav check`

PAV can be used to check the annotations of existing plastid genome GenBank files. To use this feature, provide a folder of GenBank files as input to the `pav check` command. 

### Important considerations when using existing GenBank files

When using `pav check` with existing GenBank files, be aware of potential issues related to gene annotations that span the ends of linear sequences:

#### 1. **Cross-boundary gene annotations**
Some annotation programs (including Chloë) may annotate genes that cross the ends of a linear FASTA sequence. For example:
- The 5' end of `rps19` might be annotated at the 3' end of the input sequence
- The 3' end of `rps19` might be annotated at the 5' end of the input sequence

#### 2. **Incomplete genome assemblies**
If you provide an incomplete genome assembly (e.g., a sequence that cannot be circularised due to partial assembly of a genome that is circular *in vivo*), and this sequence has:
- A partial `rps19` gene at its 3' end
- Intergenic sequence at its 5' end

This can result in a partially incorrect `rps19` annotation. This occurs because Chloë annotates CDS regions by finding corresponding ORFs. If the correct partial `rps19` ORF at the 3' end can be extended into the 5' end (intergenic sequence) before reaching an in-frame stop codon, the resulting ORF (containing intergenic sequence) will be returned as the `rps19` CDS.

#### 3. **Sequence truncation issues**
Even if the `rps19` annotation is correct (i.e., the input sequence can be circularised but happens to be linearised halfway through the `rps19` CDS), the annotation interval extends beyond the sequence length. For example:
- `rps19` is annotated with interval `156176..156472`
- The input sequence is only 156,374 bp long

This creates an issue with PAV's current implementation, which uses the Biopython SeqIO GenBank parser to extract CDS sequences. The parser cannot handle intervals that extend beyond the sequence boundaries and instead truncates the sequence to the available length (e.g., `156176..156374`).

In these scenarios, PAV will produce gene length warnings in the report files, such as:

```
['Gene rps19 length mismatch: Expected 297 != Extracted 199 for copy 1 - check this gene manually!', 'CDS length 199 is not a multiple of 3 for copy 1', "Stop codon is TGG, expected one of ['TAA', 'TAG', 'TGA'] for copy 1"]
```

### How PAV handles these issues

When annotating FASTA sequences using `pav annotate_and_check`, PAV uses a two-pass approach to avoid these issues:

1. **First annotation**: The input sequence is annotated once with Chloë
2. **Linearization**: The location of a specified gene is recovered (`psbA` by default, can be changed with `--linearise_gene <gene_name>`) and the input FASTA is linearised just upstream
3. **Second annotation**: The input sequence is annotated again with Chloë using the linearised sequence as input

This approach ensures that no annotation should cross the ends of the annotated FASTA sequence.

**Note**: If the specified gene for linearisation cannot be found, the FASTA sequence is processed as-is. For fragmented assemblies with multiple FASTA records, `psbA` may be found (and used for linearisation) in one sequence, but others will be processed without linearisation, potentially leading to the issues described above.


## Gene types supported

PAV processes and validates three main types of plastid genes:

- **CDS (Coding Sequences)**: Protein-coding genes with translation validation and backtranslation alignment
- **rRNA (Ribosomal RNA)**: Ribosomal RNA genes with direct nucleotide alignment
- **tRNA (Transfer RNA)**: Transfer RNA genes with direct nucleotide alignment

### Reference sequence naming
Reference sequences are named using the format: `{Order}_{Family}_{Genus}_{Species}_{GeneName}_{Filename}`
- Multiple copies of the same gene from a single reference file are labeled with `_copy_1`, `_copy_2`, etc.
- This ensures unique identification of each gene copy in alignment files

## File formats

### Input
- **FASTA**: Genome assemblies (.fasta, .fa, .fas) - for `pav annotate_and_check`
- **GenBank**: Annotated genomes (.gb, .gbk, .gb.gz, .gbk.gz) - for `pav check`

### Output
- **GenBank**: Annotated genomes in GenBank format (from Chloë)
- **GFF**: Gene feature format files (from Chloë)
- **EMBL**: European Molecular Biology Laboratory format
- **FASTA**: Aligned sequences and linearised genomes

## Reference data

PAV supports multiple sources of reference genomes:

### Built-in references
- **Order-specific references**: Located in [`data/order_genomes`](https://github.com/chrisjackson-pellicle/PAV/tree/main/plastid_annotation_validator/data/order_genomes) directory, organized by taxonomic order
- **Default references**: Located in [`data/reference_genomes_default`](https://github.com/chrisjackson-pellicle/PAV/tree/main/plastid_annotation_validator/data/reference_genomes_default) directory

### Custom references
- **Custom reference folder**: Users can provide their own folder of reference GenBank files using `--custom_refs_folder`
- Supports both compressed (`.gz`) and uncompressed GenBank files
- Can be used alone or in combination with built-in references

### Reference combination
PAV can combine reference sequences from multiple sources:
- Order-specific references (via `--refs_order`)
- Custom references (via `--custom_refs_folder`)
- Default references (when no specific references are specified)

All reference sources are merged to provide comprehensive validation and alignment data.

### Available orders
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

### Median lengths
Reference gene lengths are stored in [`data/plDNA_genes_median_lengths.csv`](https://github.com/chrisjackson-pellicle/PAV/blob/main/plastid_annotation_validator/data/plDNA_genes_median_lengths.csv) for validation.

## Troubleshooting

### Common issues

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
