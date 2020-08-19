#!/usr/bin/env python

'''
inStrain - parse command-line arguemnts
'''

__author__ = "Matt Olm and Alex Crits-Christoph"
__license__ = "MIT"
__email__ = "mattolm@gmail.com"
__status__ = "Development"

import os
import sys
import argparse

# Get the version
from ._version import __version__

"""
########################################
#    Argument Parsing                  #
########################################
"""

class SmartFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

def printHelp():
    print('')
    print('                ...::: inStrain v' + __version__ + ' :::...''')
    print('')
    print('  Matt Olm and Alex Crits-Christoph. MIT License. Banfield Lab, UC Berkeley. 2019''')
    print('''\

  Choose one of the operations below for more detailed help. See https://instrain.readthedocs.io for documentation.
  Example: inStrain profile -h

  Workflows:
    profile         -> Create an inStrain profile (microdiversity analysis) from a mapping.
    compare         -> Compare multiple inStrain profiles (popANI, coverage_overlap, etc.)

  Single operations:
    profile_genes   -> Calculate gene-level metrics on an inStrain profile [DEPRECATED; USE profile INSTEAD]
    genome_wide     -> Calculate genome-level metrics on an inStrain profile
    quick_profile   -> Quickly calculate coverage and breadth of a mapping using coverM
    filter_reads    -> Commands related to filtering reads from .bam files
    plot            -> Make figures from the results of "profile" or "compare"
    other           -> Other miscellaneous operations
    ''')

def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=SmartFormatter)
    subparsers = parser.add_subparsers(help='Desired operation',dest='operation')

    # Make a parent parser for all of the subparsers
    parent_parser = argparse.ArgumentParser(add_help=False)

    Bflags = parent_parser.add_argument_group('SYSTEM PARAMETERS')
    Bflags.add_argument('-p','--processes',help='Number of processes to use',default=6,type=int)
    Bflags.add_argument('-d','--debug',help='Make extra debugging output',default=False,
                        action= "store_true")
    Bflags.add_argument("-h", "--help", action="help", help="show this help message and exit")
    Bflags.add_argument(
        "--version",
        action="version",
        version="inStrain version {version}".format(version=__version__))


    # Make a parent parser for read filtering
    readfilter_parent = argparse.ArgumentParser(add_help=False)
    fiflags = readfilter_parent.add_argument_group('READ FILTERING OPTIONS')
    fiflags.add_argument("-l", "--min_read_ani", action="store", default=0.95, type=float, \
        help='Minimum percent identity of read pairs to consensus to use the reads. Must be >, not >=')
    fiflags.add_argument("--min_mapq", action="store", default=-1, type=int,\
        help='Minimum mapq score of EITHER read in a pair to use that pair. Must be >, not >=')
    fiflags.add_argument("--max_insert_relative", action="store", default=3, type=float, \
        help='Multiplier to determine maximum insert size between two reads - default is to use 3x median insert size. Must be >, not >=')
    fiflags.add_argument("--min_insert", action="store", default=50, type=int,\
        help='Minimum insert size between two reads - default is 50 bp. If two reads are 50bp each and overlap completely, their insert will be 50. Must be >, not >=')
    fiflags.add_argument("--pairing_filter", help="R|How should paired reads be handled?\n" \
        + "paired_only = Only paired reads are retained\n" \
        + 'non_discordant = Keep all paired reads and singleton reads that map to a single scaffold\n' \
        + "all_reads = Keep all reads regardless of pairing status (NOT RECOMMENDED; See documentation for deatils)\n", \
        default = "paired_only", choices={'paired_only', 'non_discordant', 'all_reads'})
    fiflags.add_argument("--priority_reads", help='The location of a list ' \
        + "of reads that should be retained regardless of pairing status " \
        + "(for example long reads or merged reads). This can be a .fastq " \
        + "file or text file with list of read names (will assume file is " \
        + "compressed if ends in .gz", default=None)

    # Make a parent parser for read output
    readoutput_parent = argparse.ArgumentParser(add_help=False)
    fiflags = readoutput_parent.add_argument_group('READ OUTPUT OPTIONS')
    # fiflags.add_argument("-s", "--generate_sam", action="store", default=None, \
    #     help='Specify the location to write a .sam file with filtered reads only.')
    fiflags.add_argument("--detailed_mapping_info", action="store_true", default=False, help='Make a detailed read report indicating deatils about each individual mapped read')

    # Make a parent parser for SNV calling
    variant_parent = argparse.ArgumentParser(add_help=False)
    fiflags = variant_parent.add_argument_group('VARIANT CALLING OPTIONS')
    fiflags.add_argument("-c", "--min_cov", action="store", default=5, type=int, \
        help='Minimum coverage to call an variant')
    fiflags.add_argument("-f", "--min_freq", action="store", default=0.05, type=float, \
        help='Minimum SNP frequency to confirm a SNV (both this AND the  FDR snp count cutoff must be true to call a SNP).')
    fiflags.add_argument("-fdr", "--fdr", action="store", default=1e-6, type=float,\
        help='SNP false discovery rate- based on simulation data with a 0.1 percent error rate (Q30)')

    # Make a parent for profile_genes
    genes_parent = argparse.ArgumentParser(add_help=False)
    Rflags = genes_parent.add_argument_group('GENE PROFILING OPTIONS')
    Rflags.add_argument("-g", "--gene_file", action="store", default=None, \
        help='Path to prodigal .fna genes file. If file ends in .gb or .gbk, will treat as a genbank file (EXPERIMENTAL; the name of the gene must be in the gene qualifier)')

    # Make a parent for genome_wide
    geneomewide_parent = argparse.ArgumentParser(add_help=False)
    Rflags = geneomewide_parent.add_argument_group('GENOME WIDE OPTIONS')
    Rflags.add_argument('-s', '--stb', help="Scaffold to bin. This can be a file with each line listing a scaffold and a bin name, tab-seperated. This can also be a space-seperated list of .fasta files, with one genome per .fasta file. If nothing is provided, all scaffolds will be treated as belonging to the same genome",
                        nargs='*', default=[])


    # Make a parent for handling mm
    mm_parent = argparse.ArgumentParser(add_help=False)
    Rflags = mm_parent.add_argument_group('READ ANI OPTIONS')
    Rflags.add_argument('--mm_level', help="Create output files on the mm level (see documentation for info)",
                        action='store_true', default=False)
    Rflags.add_argument('--skip_mm_profiling', action='store_true', default=False,\
        help="Dont perform analysis on an mm level; saves RAM and time; impacts plots and raw_data")

    '''
    ####### Arguments for profile operation ######
    '''
    # Make a parent for profile to go above the system arguments
    profile_parent = argparse.ArgumentParser(add_help=False)
    # Required positional arguments
    Rflags = profile_parent.add_argument_group('REQUIRED')
    Rflags.add_argument("bam", help="Sorted .bam file")
    Rflags.add_argument("fasta", help="Fasta file the bam is mapped to")
    # I/O Parameters
    Iflags = profile_parent.add_argument_group('I/O PARAMETERS')
    Iflags.add_argument("-o", "--output", action="store", default='inStrain', \
        help='Output prefix')
    Iflags.add_argument('--use_full_fasta_header', action='store_true', default=False,
        help='Instead of using the fasta ID (space in header before space), use the full header. Needed for some mapping tools (including bbMap)')

    profile_parser = subparsers.add_parser("profile",formatter_class=SmartFormatter,\
                    parents = [profile_parent, parent_parser, readfilter_parent, readoutput_parent, variant_parent, genes_parent, geneomewide_parent, mm_parent], add_help=False)

    # Other Parameters
    Oflags = profile_parser.add_argument_group('PROFILE OPTIONS')
    Oflags.add_argument('--database_mode', action='store_true', default=False,\
        help="Set a number of parameters to values appropriate for mapping to a " \
        + "large fasta file. Will set: --min_read_ani 0.92 --skip_mm_profiling --min_genome_coverage 1")
    Oflags.add_argument("--min_scaffold_reads", action="store", default=1, type=int,\
        help='Minimum number of reads mapping to a scaffold to proceed with profiling it')
    Oflags.add_argument("--min_genome_coverage", action="store", default=0, type=float,\
        help='Minimum number of reads mapping to a genome to proceed with profiling it. MUST profile .stb if this is set')
    Oflags.add_argument("--min_snp", action="store", default=20, \
        help='Absolute minimum number of reads connecting two SNPs to calculate LD between them.')
    Oflags.add_argument('--store_everything', action='store_true', default=False,\
        help="Store intermediate dictionaries in the pickle file; will result in significantly more RAM and disk usage")
    Oflags.add_argument("--scaffolds_to_profile", action="store",\
        help='Path to a file containing a list of scaffolds to profile- if provided will ONLY profile those scaffolds')
    Oflags.add_argument("--rarefied_coverage", action='store', default=50,\
        help='When calculating nucleotide diversity, also calculate a rarefied version with this much coverage')
    Oflags.add_argument('--window_length', action='store', default=10000, type=int,\
        help='Break scaffolds into windows of this length when profiling')

    # Other Parameters
    Iflags = profile_parser.add_argument_group('OTHER  OPTIONS')
    Iflags.add_argument('--skip_genome_wide', action='store_true', default=False,\
        help="Do not generate tables that consider groups of scaffolds belonging to genomes")
    Iflags.add_argument('--skip_plot_generation', action='store_true', default=False,\
        help="Do not make plots")


    '''
    ####### Arguments for compare operation ######
    '''
    # Make a parent for profile to go above the system arguments
    compare_parent = argparse.ArgumentParser(add_help=False)
    # Required positional arguments
    Rflags = compare_parent.add_argument_group('REQUIRED')
    Rflags.add_argument('-i', '--input', help="A list of inStrain objects, all mapped to the same .fasta file",
                        nargs='*', required=True)
    Rflags.add_argument("-o", "--output", action="store", default='instrainComparer', \
                        help='Output prefix')


    compare_parser = subparsers.add_parser("compare",formatter_class=SmartFormatter,\
                    parents = [compare_parent, parent_parser, variant_parent], add_help=False)

    # Other Parameters
    Oflags = compare_parser.add_argument_group('OTHER OPTIONS')
    Oflags.add_argument("-s", "--scaffolds", action="store", \
                        help='Location to a list of scaffolds to compare. You can also make this a .fasta file and it will load the scaffold names')
    Oflags.add_argument('--store_coverage_overlap', action='store_true', default=False,\
        help="Also store coverage overlap on an mm level")
    Oflags.add_argument('--store_mismatch_locations', action='store_true', default=False,\
        help="Store the locations of SNPs")
    Oflags.add_argument('--include_self_comparisons', action='store_true', default=False,\
        help="Also compare IS profiles against themself")

    Gflags = compare_parser.add_argument_group('GREEDY CLUSTERING OPTIONS [THIS SECTION IS EXPERIMENTAL!]')
    Gflags.add_argument('--greedy_clustering', action='store_true', default=False,\
        help="Dont do pair-wise comparisons, do greedy clustering to only find the number of clsuters. If this is set, use the parameters below as well")
    Gflags.add_argument('--g_ani', action='store', default=0.99, type=float,\
        help="ANI threshold for greedy clustering- put the fraction not the percentage (e.g. 0.99, not 99)")
    Gflags.add_argument('--g_cov', action='store', default=0.99, type=float,\
        help="Alignment coverage for greedy clustering- put the fraction not the percentage (e.g. 0.5, not 10)")
    Gflags.add_argument('--g_mm', action='store', default=100, type=int,\
        help="Maximum read mismatch level")

    '''
    ####### Arguments for profile_genes operation ######
    '''
    # Make a parent for profile to go above the system arguments
    genes_io = argparse.ArgumentParser(add_help=False)

    # Required positional arguments
    Rflags = genes_io.add_argument_group('INPUT / OUTPUT')
    Rflags.add_argument("-i", '--IS', help="an inStrain profile object", required=True)
    Rflags.add_argument('--store_everything', action='store_true', default=False,\
        help="Store gene sequences in the IS object")

    genes_parser = subparsers.add_parser("profile_genes",formatter_class=SmartFormatter,\
                    parents = [genes_parent, genes_io, parent_parser], add_help=False)

    '''
    ####### Arguments for genome_wide operation ######
    '''
    # Make a parent for profile to go above the system arguments
    genome_parser = subparsers.add_parser("genome_wide",formatter_class=SmartFormatter,\
                    parents = [geneomewide_parent, genes_io, mm_parent, parent_parser], add_help=False)

    '''
    ####### Arguments for plot operation ######
    '''
    # Make a parent for profile to go above the system arguments
    plot_parent = argparse.ArgumentParser(add_help=False)

    # Required positional arguments
    Rflags = plot_parent.add_argument_group('REQUIRED')
    Rflags.add_argument("-i", '--IS', help="an inStrain profile object", required=True)
    Rflags.add_argument("-pl", "--plots", help= "R|Plots. "
                        + "Input 'all' or 'a' to plot all\n"
                        + "1) Coverage and breadth vs. read mismatches\n"
                        + "2) Genome-wide microdiversity metrics\n"
                        + "3) Read-level ANI distribution\n"
                        + "4) Major allele frequencies\n"
                        + "5) Linkage decay\n"
                        + "6) Read filtering plots\n"
                        + "7) Scaffold inspection plot (large)\n"
                        + "8) Linkage with SNP type (GENES REQUIRED)\n"
                        + "9) Gene histograms (GENES REQUIRED)\n"
                        + "10) Compare dendrograms (RUN ON COMPARE; NOT PROFILE)\n",
                        nargs='*', default='a')

    POflags = plot_parent.add_argument_group('OPTIONAL FIGURE ADJUSTMENTS')
    POflags.add_argument("-mb", "--minimum_breadth", default=0, type=float,
                    help= "Minimum breadth of coverage for genome to make it into plot (from 0-1).")
    POflags.add_argument("-g", "--genomes", nargs='*',
                    help= "Only plot genomes with the names provided in this argument")

    plot_parser = subparsers.add_parser("plot",formatter_class=SmartFormatter,\
                    parents = [plot_parent, parent_parser], add_help=False)

    '''
    ####### Arguments for quick_profile operation ######
    '''
    # Make a parent for profile to go above the system arguments
    quick_parent = argparse.ArgumentParser(add_help=False)
    # Required positional arguments
    Rflags = quick_parent.add_argument_group('REQUIRED')
    Rflags.add_argument("bam", help="Sorted .bam file")
    Rflags.add_argument("fasta", help="Fasta file the bam is mapped to")


    quick_parser = subparsers.add_parser("quick_profile",formatter_class=SmartFormatter,\
                    parents = [quick_parent, parent_parser], add_help=False)

    # Other Parameters
    Oflags = quick_parser.add_argument_group('OTHER OPTIONS')
    Oflags.add_argument('-s', '--stb', help="Scaffold to bin. This can be a file with each line listing a scaffold and a bin name, tab-seperated. This can also be a space-seperated list of .fasta files, with one genome per .fasta file. If nothing is provided, all scaffolds will be treated as belonging to the same genome",
                        nargs='*', default=[])
    Oflags.add_argument("-o", "--output", action="store", \
                        help='Output prefix', default='QuickProfile')
    Oflags.add_argument("--breadth_cutoff", type=float, default=0.5,
                        help='Minimum genome breadth to pull scaffolds')
    Oflags.add_argument("--stringent_breadth_cutoff", type=float, default=0.00,
                        help='Minimum breadth to let scaffold into coverm raw results (done with greater than; NOT greater than or equal to)')


    '''
    ####### Arguments for filter_reads operation ######
    '''
    # Make a parent for profile to go above the system arguments
    reads_parent = argparse.ArgumentParser(add_help=False)
    # Required positional arguments
    Rflags = reads_parent.add_argument_group('REQUIRED')
    Rflags.add_argument("bam", help="Sorted .bam file")
    Rflags.add_argument("fasta", help="Fasta file the bam is mapped to")
    Rflags.add_argument("-o", "--output", action="store", \
                        help='Location of folder to store read report(s)')


    reads_parser = subparsers.add_parser("filter_reads",formatter_class=SmartFormatter,\
                    parents = [reads_parent, parent_parser, readfilter_parent,
                    readoutput_parent], add_help=False)

    '''
    ####### Arguments for other operation ######
    '''
    other_parser = subparsers.add_parser("other",formatter_class=SmartFormatter,\
                    parents = [parent_parser], add_help=False)

    # Other Parameters
    Oflags = other_parser.add_argument_group('OTHER OPTIONS')
    Oflags.add_argument('--old_IS', help="Convert an old inStrain version object to the newer version.")
    Oflags.add_argument('--run_statistics', help='Generate runtime reports for an inStrain run.')

    '''
    ####### PARSE THE ARGUMENTS ######
    '''

    # Handle the situation where the user wants the raw help
    if (len(args) == 0 or args[0] == '-h' or args[0] == '--help'):
        printHelp()
        sys.exit(0)
    else:
        return parser.parse_args(args)
