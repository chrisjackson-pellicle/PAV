#!/usr/bin/env python

"""
Contains subparsers used in pav.py
"""


def _add_gene_length_warnings_group(parser):
    """Add gene length warnings argument group to parser"""
    optional_group_gene_length_warnings = parser.add_argument_group('Gene length warnings')

    optional_group_gene_length_warnings.add_argument('--min_length_percentage', '-min_len',
                                type=float,
                                default=0.8,
                                metavar='FLOAT',
                                dest='min_length_percentage',
                                help='Minimum length percentage of the gene median length for a warning to be issued. Default is: %(default)s')
    optional_group_gene_length_warnings.add_argument('--max_length_percentage', '-max_len',
                                type=float,
                                default=1.2,
                                metavar='FLOAT',
                                dest='max_length_percentage',
                                help='Maximum length percentage of the gene median length for a warning to be issued. Default is: %(default)s')


def _add_alignment_group(parser):
    """Add alignment with reference genes argument group to parser"""
    optional_group_alignment = parser.add_argument_group('Alignment with reference genes')

    optional_group_alignment.add_argument('--no_alignment', '-no_align',
                                action='store_true',
                                default=False,
                                dest='no_alignment',
                                help='Do not align annotated genes with reference genes. Default is: %(default)s')
    optional_group_alignment.add_argument('--refs_order', '-refs_ord',
                                metavar='ORDER',
                                type=str,
                                nargs='+',
                                default=[],
                                dest='refs_order',
                                help='Order(s) to use for reference genes. Can be specified multiple times. Default is: %(default)s')
    optional_group_alignment.add_argument('--custom_refs_folder', '-custom_refs',
                                metavar='DIR',
                                type=str,
                                default=None,
                                dest='custom_refs_folder',
                                help='Custom folder containing reference GenBank files. Can be used in addition to --refs_order or default references.')


def _add_intergenic_group(parser):
    """Add intergenic region analysis argument group to parser"""
    intergenic_group = parser.add_argument_group('Intergenic region analysis')
    
    intergenic_group.add_argument('--min_intergenic_length', '-min_ig_len',
                                 metavar='INTEGER',
                                 type=int,
                                 default=0,
                                 dest='min_intergenic_length',
                                 help='Minimum length of intergenic region to analyze. Default is: %(default)s')
    
    intergenic_group.add_argument('--blast_evalue', '-evalue',
                                 metavar='FLOAT',
                                 type=float,
                                 default=1e-10,
                                 dest='blast_evalue',
                                 help='BLAST E-value threshold for intergenic region analysis. Default is: %(default)s')
    
    intergenic_group.add_argument('--skip_intergenic_analysis', '-skip_ig',
                                 action='store_true',
                                 default=False,
                                 dest='skip_intergenic_analysis',
                                 help='Skip intergenic region analysis. Default is: %(default)s')
    
    intergenic_group.add_argument('--debug_intergenic', '-debug_ig',
                                 action='store_true',
                                 default=False,
                                 dest='debug_intergenic',
                                 help='Write intergenic regions to FASTA files for debugging. Default is: %(default)s')
    
    intergenic_group.add_argument('--max_blast_hits', '-max_hits',
                                 metavar='INTEGER',
                                 type=int,
                                 default=1,
                                 dest='max_blast_hits',
                                 help='Maximum number of BLAST hits to retain per intergenic region. Default is: %(default)s')
    
    intergenic_group.add_argument('--custom_blast_db', '-blast_db',
                                 metavar='PATH',
                                 type=str,
                                 default=None,
                                 dest='custom_blast_db',
                                 help='Custom BLAST database path or directory for intergenic region analysis. Can be a direct path to the database or a directory containing BLAST database files. If not provided, uses the default order_genomes_blastdb.')


def _add_general_options_group(parser, include_chloe=True):
    """Add general pipeline options argument group to parser"""
    optional_group_general = parser.add_argument_group('General pipeline options')
    optional_group_general.add_argument('--output_directory', '-out_dir',
                                type=str,
                                default='output_directory',
                                metavar='DIR',
                                dest='output_directory',
                                help='Output directory for the subcommand. Default is: %(default)s')
    optional_group_general.add_argument('--pool', '-p',
                                 type=int,
                                 default=1,
                                 metavar='INTEGER',
                                 dest='pool',
                                 help='The number of CPUs to use for the subcommand. Default is: %(default)s')
    optional_group_general.add_argument('--threads', '-t',
                                 type=int,
                                 default=1,
                                 metavar='INTEGER',
                                 dest='threads',
                                 help='The number of threads to use for the subcommand. Default is: %(default)s')
    
    optional_group_general.add_argument('--run_profiler',
                              action='store_true',
                              dest='run_profiler',
                              default=False,
                              help='If supplied, run the subcommand using cProfile. Saves a *.csv file of results. '
                                   'Default is: %(default)s')
    
    if include_chloe:

        optional_group_general.add_argument('--linearise_gene',
                                 type=str,
                                 nargs='+',
                                 default=['psbA'],
                                 metavar='GENE_NAME',
                                 help='Gene(s) to use for genome linearisation. Can specify multiple genes. '
                                      'Genes will be tried in order until one is found. Default is: %(default)s')
    

def annotate_and_check_parser(subparsers):
    """
    Parser for annotate_and_check.py

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_annotate_and_check = subparsers.add_parser('annotate_and_check', help='Annotate and check plastid DNA')
    
    ####################################################################################################################
    required_positional_group = parser_annotate_and_check.add_argument_group('Required input')
    
    required_positional_group.add_argument('genome_fasta_dir',
                                metavar='GENOME_FASTA_DIR',
                                help='Directory containing plastid DNA FASTA files.')
    
    required_positional_group.add_argument('metadata_tsv',
                                metavar='METADATA_TSV',
                                help='TSV file containing sample metadata for EMBL conversion. Required file should contain columns: '
                                'input_filename, project_id, locus_tag, genus_species, linear_or_circular. '
                                'ALL samples must be listed. Only input_filename and linear_or_circular require values; '
                                'empty optional fields will use defaults. linear_or_circular must be "linear" or "circular".')

    required_positional_group.add_argument('chloe_project_dir',
                                metavar='CHLOE_PROJECT_DIR',
                                help='Path to the Chloë project directory; chloe.jl is expected at <DIR>/chloe.jl')

    # Add all argument groups using helper functions
    _add_gene_length_warnings_group(parser_annotate_and_check)
    _add_alignment_group(parser_annotate_and_check)
    _add_intergenic_group(parser_annotate_and_check)
    _add_general_options_group(parser_annotate_and_check, include_chloe=True)


    return parser_annotate_and_check


def check_parser(subparsers):
    """
    Parser for check.py (continues pipeline from annotated GenBank files)

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_check = subparsers.add_parser('check', help='Continue pipeline from annotated GenBank files')
    
    ####################################################################################################################
    required_positional_group = parser_check.add_argument_group('Required input')
    
    required_positional_group.add_argument('annotated_genbank_dir',
                                metavar='ANNOTATED_GENBANK_DIR',
                                help='Directory containing already annotated GenBank files.')
    
    # Add optional arguments
    optional_group_general = parser_check.add_argument_group('Optional input')
    optional_group_general.add_argument('--metadata_tsv',
                                metavar='METADATA_TSV',
                                default=None,
                                help='TSV file containing sample metadata for EMBL conversion. If not provided, samples '
                                'will be processed without metadata.')


    # Add all argument groups using helper functions (excluding chloë options)
    _add_gene_length_warnings_group(parser_check)
    _add_alignment_group(parser_check)
    _add_intergenic_group(parser_check)
    _add_general_options_group(parser_check, include_chloe=False)

    return parser_check


def db_for_intgen_parser(subparsers):
    """
    Parser for db_for_integen.py

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_db_for_intgen = subparsers.add_parser('db_for_intgen', help='Extract genes and features from GenBank files and create BLAST database for intergenic analyses')
    
    ####################################################################################################################
    required_positional_group = parser_db_for_intgen.add_argument_group('Required input')
    
    required_positional_group.add_argument('input_dir',
                                metavar='INPUT_DIR',
                                help='Input directory containing GenBank files (.gb, .gbk, .gb.gz, or .gbk.gz)')
    
    ####################################################################################################################
    optional_group = parser_db_for_intgen.add_argument_group('Optional parameters')
    
    optional_group.add_argument('--output_directory', '-out_dir',
                                type=str,
                                default='output_directory_intgen_db',
                                dest='output_directory',
                                metavar='DIR',
                                help='Output directory for the subcommand. Default is: %(default)s')
    optional_group.add_argument('--min_length', 
                               type=int, 
                               default=10,
                               metavar='INTEGER',
                               help='Minimum sequence length to include. Default is: %(default)s')
    optional_group.add_argument('--run_profiler',
                              action='store_true',
                              dest='run_profiler',
                              default=False,
                              help='If supplied, run the subcommand using cProfile. Saves a *.csv file of results. '
                                   'Default is: %(default)s')

    return parser_db_for_intgen