### Version 1.0.1
- Added input FASTA validation to ensure all sequences contain only valid DNA characters (including ambiguity codes, Ns, dashes, and whitespace)
- Fixed GenBank topology issue: Chloë always outputs "circular" in GenBank LOCUS lines, but PAV now corrects this to match the topology specified in the metadata file (`linear_or_circular` field) when running `pav annotate_and_check`
- Restructured output directories. Renamed:
  - `03_alignments_with_refs` -> `03_alignments`. Subdirectories renamed: 
  - `01_per_sample_alignments` -> `01_per_sample_alignments_with_refs`
  - `02_per_gene_alignments` -> `02_per_gene_alignments_with_refs`
- PAV now generates sample-only alignments (no reference sequences) in addition to standard reference-based alignments. These are written to `03_alignments/03_per_gene_alignments_sample_only`
- Removed `--no_alignment` parameter
- Changed alignment behavior: all alignment types now run by default with new options to disable specific types:
  - `--no_per_sample_alignments_with_refs`: Disables per-sample alignments with references
  - `--no_per_gene_alignments_with_refs`: Disables per-gene alignments with references (combined samples)
  - `--no_sample_only_alignments`: Disables sample-only alignments
- Improved error handling when running Chloë so that errors for a single sample/sequence are logged and do not crash the whole pipeline
- Added `--no-filter` option to `pav annotate_and_check`to disable Chloë annotation filtering. When present, the `--no-filter` flag will be used in the Chloë annotation command. Previously this flag was used by default
- Fixed debug intergenic FASTA file output: when using `--debug_intergenic` with `pav check`, debug FASTA files are now written to the `04_intergenic_analysis` output directory instead of sometimes overwriting input GenBank files. These FASTA files are also now written to `04_intergenic_analysis` when running `pav annotate_and_check`
- Added Chloë installation verification: PAV now checks if Chloë is properly installed and accessible by running `chloe.jl annotate --help` before starting the annotation pipeline
- Enhanced error handling: PAV now exits gracefully with an informative message if no annotated genomes are found after Chloë processing, preventing downstream errors
- Added cleaned/adjusted GenBank file output: PAV now generates cleaned GenBank files, concatentated for multi-sequence samples. These files are produced using the same cleaning logic as EMBL conversion (locus tags, source features, etc.)
- Added final per-sample files output for Genbank/GFF/FASTA: PAV now generates organised final files in `01_annotated_genomes/01_final_per_sample_files/` directory with per-sample subdirectories containing:
  - `{sample_name}_final.gbk`: Cleaned/adjusted GenBank files, concatenated for multi-sequence samples
  - `{sample_name}_final.gff`: Concatenated GFF files (direct output of Chloë for the moment)
  - `{sample_name}_final.fasta`: Final FASTA files (linearised version if produced, concatenated for multi-sequence samples)


### Version 1.0.0
- Optional parameters `--chloe_project_dir` and `--chloe_project_name` have been removed. PAV now expects the path to the Chloe project directory as a required positional argument to `pav annotate_and_check`. PAV no longer checks for the Chloe project directory in $CONDA_PREFIX/bin
- Added path resolution for the PAV data directory, as setup for a conda package install
- Set 'spawn' as multiprocessing start method for Linux, to avoid listener thread locks causing program hang
- Capture KeyboardInterrupt when running program to ensure listener thread shutdown on exit

### Version 0.0.1 (Initial Release)
- **Initial release** with plastid genome annotation and validation functions
- **Genome annotation** with Chloë integration
- **Reference-based validation** against multiple taxonomic orders
- **Alignment generation** for CDS, rRNA, and tRNA genes
- **Quality assessment** with gene length checking and translation validation
- **EMBL/ENA template generation** for sequence submission
- **Intergenic region analysis** with BLAST database support
- **Multi-sequence FASTA support** for fragmented assemblies
- **Parallel processing** for improved performance
- **Comprehensive reporting** and logging system