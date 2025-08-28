"""
Plastid Annotation Validator (PAV) - Main Entry Point

PAV is a comprehensive tool for annotating and validating plastid genome annotations.
This module serves as the main command-line interface and entry point for the PAV pipeline.

Overview:
    PAV automates the process of annotating plastid genomes using the chloë annotation
    pipeline and provides extensive validation capabilities including gene length checking,
    multiple sequence alignment generation, and intergenic region analysis.

Main Features:
    - Automated plastid genome annotation using chloë
    - Gene annotation validation against reference data
    - Multiple sequence alignment generation (CDS, rRNA, per-gene)
    - Intergenic region extraction and BLAST analysis
    - Comprehensive reporting and logging
    - Parallel processing for improved performance
    - GenBank feature extraction and BLAST database creation for intergenic region analysis

Usage:
    pav annotate_and_check [options]
    pav check [options]
    pav db_for_intgen [options]
    
    For help on specific subcommands:
    pav annotate_and_check --help
    pav check --help
    pav db_for_intgen --help

Subcommands:
    annotate_and_check    Main annotation and validation pipeline
    check                 Continue pipeline from annotated GenBank files
    db_for_intgen         Extract genes and features from GenBank files and create BLAST database for intergenic region analysis
    
Options:
    --version            Show version information
    --help               Show help message
    --run_profiler       Enable performance profiling

Examples:
    # Basic annotation and validation
    pav annotate_and_check genomes/ metadata.tsv -o output/
    
    # With custom parameters
    pav annotate_and_check genomes/ metadata.tsv -o output/ --min-length-percentage 80 --max-length-percentage 120
    
    # Continue pipeline from annotated GenBank files
    pav check annotated_genbank/ metadata.tsv -o output/
    
    # Extract features from GenBank files and create BLAST database for intergenic region analysis
    pav db_for_intgen genbank_files/ -o output/ --min_length 50

Dependencies:
    - chloë: Genome annotation pipeline
    - MAFFT: Multiple sequence alignment
    - BLAST+: Sequence similarity search
    - TrimAl: Alignment trimming
    - Biopython: Sequence manipulation

Author: Chris Jackson chris.jackson@rbg.vic.gov.au
License: See LICENSE file
"""

import textwrap
import sys
import cProfile
import argparse
import os

# Import plastid_annotation_validator modules:
import plastid_annotation_validator.utils as utils
import plastid_annotation_validator.annotate_and_check as annotate_and_check
import plastid_annotation_validator.pav_subparsers as pav_subparsers
import plastid_annotation_validator.db_for_intgen as db_for_intgen
from plastid_annotation_validator.version import __version__


def annotate_and_check_main(args):
    """
    Calls the function main() from module annotate_and_check

    :param args: argparse namespace with subparser options for function annotate_and_check()
    :return: None: no return value specified; default is None
    """

    if len(sys.argv) > 1:
        subcommand_name = sys.argv[1]
    else:
        subcommand_name = 'annotate_and_check'

    args.subcommand_name = subcommand_name

    # Run main() from module annotate_and_check.py
    annotate_and_check.main(args)


def check_main(args):
    """
    Calls the function check_pipeline() from module annotate_and_check

    :param args: argparse namespace with subparser options for function check_pipeline()
    :return: None: no return value specified; default is None
    """

    if len(sys.argv) > 1:
        subcommand_name = sys.argv[1]
    else:
        subcommand_name = 'check'

    args.subcommand_name = subcommand_name

    # Run check_pipeline() from module annotate_and_check
    annotate_and_check.check_pipeline(args)


def db_for_intgen_main(args):
    """
    Calls the function main() from module db_for_intgen

    :param args: argparse namespace with subparser options for function db_for_intgen.main()
    :return: None: no return value specified; default is None
    """

    if len(sys.argv) > 1:
        subcommand_name = sys.argv[1]
    else:
        subcommand_name = 'db_for_intgen'

    args.subcommand_name = subcommand_name

    # Run main() from module db_for_intgen.py
    db_for_intgen.main(args)


def parse_arguments():
    """
    Creates main parser and add subparsers. Parses command line arguments

    :return argparse.Namespace arguments: arguments for the given command/subcommand
    """

    parser = argparse.ArgumentParser(prog='pav',
                                     description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='To view parameters and help for a subcommand, use e.g. "pav '
                                            '--help"')
    parser.add_argument('--version', '-v',
                        action='version', 
                        version=__version__,
                        help='Print the pav version number.')

    # Add subparsers:
    subparsers = parser.add_subparsers(title='Subcommands for PAV', metavar='')
    parser_annotate_and_check = pav_subparsers.annotate_and_check_parser(subparsers)
    parser_check = pav_subparsers.check_parser(subparsers)
    parser_db_for_intgen = pav_subparsers.db_for_intgen_parser(subparsers)

    # Set functions for subparsers:
    parser_annotate_and_check.set_defaults(func=annotate_and_check_main)
    parser_check.set_defaults(func=check_main)
    parser_db_for_intgen.set_defaults(func=db_for_intgen_main)

    # Parse arguments
    arguments = parser.parse_args()

    # If no subcommand was provided, show help
    if not hasattr(arguments, 'func'):
        parser.print_help()
        sys.exit(0)

    return arguments


def main():

    title = textwrap.dedent(
        fr"""
                ██████╗      █████╗     ██╗   ██╗
                ██╔══██╗    ██╔══██╗    ██║   ██║
                ██████╔╝    ███████║    ██║   ██║
                ██╔═══╝     ██╔══██║    ╚██╗ ██╔╝
                ██║         ██║  ██║     ╚████╔╝ 
                ╚═╝         ╚═╝  ╚═╝      ╚═══╝  

        """
    )

    sys.stdout.write(title)
    sys.stdout.flush()

    # Parse arguments for the command/subcommand used:
    args = parse_arguments()

    # Run the function associated with the subcommand, with or without cProfile:
    if args.run_profiler:
        profiler = cProfile.Profile()
        profiler.enable()
        args.func(args)
        profiler.disable()
        csv = utils.cprofile_to_csv(profiler)
        
        with open(f'{sys.argv[1]}_cprofile.csv', 'w+') as cprofile_handle:
            cprofile_handle.write(csv)
    else:
        
        args.func(args)

if __name__ == '__main__':
    main()