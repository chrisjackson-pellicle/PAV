### Version 1.0.0

- Optional parameters `--chloe_project_dir` and `--chloe_project_name` have been removed. PAV now expects the path to the Chloe project directory as a required positional argument to `pav annotate_and_check`. PAV no longer checks for the Chloe project directory in $CONDA_PREFIX/bin.
- Added path resolution for the PAV data directory, as setup for a conda package install.
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