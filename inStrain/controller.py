# Get the version
from ._version import __version__

# Import packages
import os
import sys
import h5py
import logging
import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# Import inStrain stuff
import inStrain.profileUtilities
import inStrain.filter_reads
import inStrain.readComparer
import inStrain.GeneProfile
import inStrain.genomeUtilities
import inStrain.plottingUtilities
import inStrain.quickProfile
import inStrain.SNVprofile

class Controller():
    '''
    Controller of the whole shebang
    '''
    def main(self, args):
        ''' Parse user options and call the correct pipeline'''
        # Call the appropriate workflow
        if args.operation == "profile":
            self.profile_operation(args)

        if args.operation == "compare":
            self.compare_operation(args)

        if args.operation == "filter_reads":
            self.filter_reads_operation(args)

        if args.operation == "profile_genes":
            self.profile_genes_operation(args)

        if args.operation == "genome_wide":
            self.genome_wide_operation(args)

        if args.operation == "quick_profile":
            self.quick_profile_operation(args)

        if args.operation == "plot":
            self.plot_operation(args)

        if args.operation == "other":
            self.other_operation(args)

    def profile_operation(self, args):
        ProfileController().main(args)

    def compare_operation(self, args):
        inStrain.readComparer.main(args)

    def filter_reads_operation(self, args):
        inStrain.filter_reads.Controller().main(args)

    def profile_genes_operation(self, args):
        inStrain.GeneProfile.Controller().main(args)

    def genome_wide_operation(self, args):
        inStrain.genomeUtilities.Controller().main(args)

    def quick_profile_operation(self, args):
        inStrain.quickProfile.main(args)

    def plot_operation(self, args):
        inStrain.plottingUtilities.main(args)

    def other_operation(self, args):
        # Check if you should convert IS profile
        if args.old_IS != None:
            inStrain.SNVprofile.convert_SNVprofile(args.old_IS)

class ProfileController():
    '''
    Main controller of the profile command
    '''

    def main(self, args):
        '''
        The main method when run on the command line
        '''

        # Parse arguments
        #args = parse_arguments(sys_args)
        if args.output == 'inStrain':
            args.output = args.fasta.split(".")[0].split("/")[-1] #default prefix is now fasta prefix -alexcc 5/8/2019

        args = self.validate_arguments(args)
        message = """\
***************************************************
    ..:: inStrain profile Step 1. Filter reads ::..
***************************************************
        """
        logging.info(message)

        global s2l # make ths global so we can access it later.

        # Parse the args
        bam = args.bam
        vargs = vars(args)
        del vargs['bam']

        # Set up .fasta file
        FAdb, s2s = self.load_fasta(args)

        # Load dictionary of paired reads
        scaffolds = list(FAdb['scaffold'].unique())
        Rdic, RR = inStrain.filter_reads.load_paired_reads2(bam, scaffolds, **vargs)
        logging.info("{0:,} read pairs remain after filtering".format(RR['filtered_pairs'].tolist()[0]))

        if RR['filtered_pairs'].tolist()[0] == 0:
            logging.error("Because no read pairs remain I'm going to crash now. Maybe this is failing because you dont have paired reads (in which case you should adjust --pairing_filter option), or maybe its failing because the mapper you used uses full fasta headers (in which case you should use the flag --use_full_fasta_header)")
            raise Exception('No paired reads detected; see above message')

        # Get scaffold to paired reads (useful for multiprocessing)
        s2p = RR.set_index('scaffold')['filtered_pairs'].to_dict()

        # Filter the .fasta file
        FAdb = self.filter_fasta(FAdb, RR, args.min_fasta_reads)
        FAdb['filtered_pairs'] = FAdb['scaffold'].map(s2p)
        FAdb.sort_values('filtered_pairs', inplace=True, ascending=False)

        if args.skip_mm_profiling:
            Rset = set(Rdic.keys())
            del Rdic
            Rdic = Rset

        if args.debug:
            for att in ['Rdic', 'RR', 'FAdb', 's2s']:
                logging.debug("RAM PROFILE: reads {0} {1:.2f} Mb".format(att,
                    sys.getsizeof(eval(att))/1e6))

        # Profile the .bam file

        vargs['s2s'] = s2s
        vargs['s2p'] = s2p

        message = """\
***************************************************
.:: inStrain profile Step 2. Profile scaffolds ::..
***************************************************
        """
        logging.info(message)


        Sprofile = inStrain.profileUtilities.profile_bam(bam, FAdb, Rdic, **vargs)

        # Add the read report
        Sprofile.store('read_report', RR, 'pandas', "Report on reads")

        if args.store_everything:
            if args.skip_mm_profiling:
                Sprofile.store('Rdic', Rdic, 'pickle', 'list of filtered read pairs')
            else:
                Sprofile.store('Rdic', Rdic, 'dictionary', 'Read pair -> mismatches')

        # Save
        logging.info('Storing output')
        self.write_output(Sprofile, args)

        # Run the rest of things
        args.IS = Sprofile.location

        # See if you can profile genes as well
        message = """\
***************************************************
  .:: inStrain profile Step 3. Profile genes ::..
***************************************************
        """
        logging.info(message)
        if args.gene_file != None:
            args.IS = Sprofile.location
            Controller().profile_genes_operation(args)
        else:
            logging.info('Nevermind! You didnt include a genes file')

        # Make things genome_wide
        message = """\
***************************************************
.:: inStrain profile Step 4. Make genome-wide ::..
***************************************************
        """
        logging.info(message)
        if not args.skip_genome_wide:
            args.IS = Sprofile.location
            Controller().genome_wide_operation(args)
        else:
            logging.info('Nevermind! You chose to skip genome_wide')

        # Generate plots
        message = """\
***************************************************
 .:: inStrain profile Step 5. Generate plots ::..
***************************************************
        """
        logging.info(message)
        if not args.skip_plot_generation:
            args.IS = Sprofile.location
            args.plots = 'a'
            Controller().plot_operation(args)
        else:
            logging.info('Nevermind! You chose to skip making plots')

        # Final message
        message = """\
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ..:: inStrain profile finished ::..

Output tables........ {0}
Figures.............. {1}

See documentation for output descriptions - https://instrain.readthedocs.io/en/latest/

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        """.format(Sprofile.get_location('output'), \
            Sprofile.get_location('figures'))
        logging.info(message)

        return Sprofile

    def load_fasta(self, args):
        '''
        Load the sequences to be profiled

        Return a table listing scaffold name, start, end
        '''
        # PROFILE ALL SCAFFOLDS IN THE .FASTA FILE
        if args.use_full_fasta_header:
            scaff2sequence = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"), key_function=_get_description)
        else:
            scaff2sequence = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta")) # set up .fasta file

        s2l = {s:len(scaff2sequence[s]) for s in list(scaff2sequence.keys())} # Get scaffold2length
        Fdb = pd.DataFrame(list(s2l.items()), columns=['scaffold', 'end'])
        Fdb['start'] = 0

        if args.scaffolds_to_profile != None:
            print(args.scaffolds_to_profile)
            Fdb = Fdb[Fdb['scaffold'].isin(args.scaffolds_to_profile)]
            s2s = {scaff:seq for scaff, seq in scaff2sequence.items() if scaff in args.scaffolds_to_profile}

        if len(Fdb) == 0:
            logging.error("The provided scaffold list has no overlap with the provided .fasta file!")
            logging.error("Example scaffolds in list: {0}".format("\n".join(args.scaffolds_to_profile)))
            sys.exit()

        return Fdb, scaff2sequence # also return s2l - alexcc 5/9/2019: Nah, make it scaff2sequence (s2s) (M.O. 6/10/19)

    def load_paired_reads(self, args, FAdb):
        '''
        Load paired reads to be profiled

        Return a dictionary of read pair -> info
        '''
        # Load paired reads
        scaffolds = list(FAdb['scaffold'].unique())
        scaff2pair2info, scaff2total = inStrain.filter_reads.get_paired_reads_multi(args.bam, scaffolds, processes=args.processes, ret_total=True)

        # Merge into single
        pair2info = {}
        for scaff, p2i in scaff2pair2info.items():
            for p, i in p2i.items():
                pair2info[p] = i

        # Make read report
        logging.info('Making read report')
        RR = inStrain.filter_reads.makeFilterReport(scaff2pair2info, scaff2total, pair2info=pair2info, **vars(args))

        # Filter the dictionary
        logging.info('Filtering reads')
        pair2infoF = inStrain.filter_reads.filter_paired_reads_dict(pair2info,
            **vars(args))

        # Make a report on these pair reads
        #make_read_report(pair2info, pair2infoF, args)

        return pair2infoF, RR

    def validate_arguments(self, args):
        '''
        Do some parsing, start up a logger
        '''
        # Set up "base"
        out_base = args.output

        # Set up Logger
        outbase = out_base
        RCprof = inStrain.SNVprofile.SNVprofile(outbase)
        log_loc = RCprof.get_location('log') + 'log.log'
        setup_logger(log_loc)

        # Make the bam file if you need to
        args.bam = inStrain.profileUtilities.prepare_bam_fie(args)

        # Load the list of scaffolds
        args.scaffolds_to_profile = load_scaff_list(args.scaffolds_to_profile)

        # Fix the fdr
        if args.fdr == 0:
            args.fdr = 1e-6

        return args

    def write_output(self, Sprofile, args):
        '''
        Write output files
        '''
        logging.debug("Writing output files now")
        out_base = Sprofile.get_location('output') + os.path.basename(Sprofile.get('location')) + '_'

        # Write the scaffold profile
        db = Sprofile.get_nonredundant_scaffold_table()
        db.to_csv(out_base + 'scaffold_info.tsv', index=False, sep='\t')

        # Write the SNP frequencies
        db = Sprofile.get_nonredundant_snv_table()
        db.to_csv(out_base + 'SNVs.tsv', index=False, sep='\t')

        # Write the linkage
        db = Sprofile.get_nonredundant_linkage_table()
        db.to_csv(out_base + 'linkage.tsv', index=False, sep='\t')

        # Write the read report
        RR = Sprofile.get('read_report')
        inStrain.filter_reads.write_read_report(RR, out_base + 'read_report.tsv', **vars(args))

    def filter_fasta(self, FAdb, RR, min_reads=0):
        '''
        Filter the .fasta file based on the min number of mapped paired reads
        '''
        s2r = RR.set_index('scaffold')['filtered_pairs'].to_dict()
        FAdb = FAdb[[True if (s2r[s] >= min_reads) else False for s in FAdb['scaffold']]]
        return FAdb

def load_scaff_list(list):
    '''
    If this is a text file of scaffolds, return it

    If it's a .fasta file, return a list of scaffolds
    '''

    if list == None:
        return None

    # Try as it its a fasta file
    scaffs = []
    handle = open(list, "r")
    fasta = SeqIO.parse(handle, "fasta")
    for f in fasta:
        scaffs.append(f.id)

    if len(scaffs) > 0:
        handle.close()
        return(set(scaffs))

    else:
        scaffs = []
        handle.close()
        handle = open(list, "r")
        for line in handle.readlines():
            scaffs.append(line.strip())
        handle.close()
        return set(scaffs)

def setup_logger(loc):
    ''' set up logger such that DEBUG goes only to file, rest go to file and console '''

    # Cancel if a logger already exists:
    if logging.getLogger('').handlers:
        return

    # set up logging everything to file
    logging.basicConfig(level=logging.DEBUG,
                       format='%(asctime)s %(levelname)-8s %(message)s',
                       datefmt='%m-%d %H:%M',
                       filename=loc)

    # set up logging of INFO or higher to sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    console.setFormatter(formatter)

    logging.getLogger('').addHandler(console)

    logging.debug("!"*80)
    logging.debug("***Logger started up at {0}***".format(loc))
    logging.debug("Command to run inStrain was: {0}\n".format(' '.join(sys.argv)))
    logging.debug("inStrain version {0} was run \n".format(__version__))
    logging.debug("!"*80 + '\n')

def _get_description(rec):
    return rec.description

# def parse_arguments(args):
#     parser = argparse.ArgumentParser(description="inStrain version {0}".format(__version__),
#              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#
#     # Required positional arguments
#     parser.add_argument("bam", help="Sorted .bam file")
#     parser.add_argument("fasta", help="Fasta file the bam is mapped to")
#
#     # Optional arguments
#     parser.add_argument("-o", "--output", action="store", default='inStrain', \
#         help='Output prefix')
#     parser.add_argument("-p", "--processes", action="store", default=6, type=int, \
#         help='Threads to use for multiprocessing')
#     parser.add_argument("-c", "--min_cov", action="store", default=5, \
#         help='Minimum SNV coverage')
#     parser.add_argument("-s", "--min_snp", action="store", default=20, \
#         help='Absolute minimum number of reads connecting two SNPs to calculate LD between them.')
#     parser.add_argument("-f", "--min_freq", action="store", default=0.05, \
#         help='Minimum SNP frequency to confirm a SNV (both this AND the 0.  001 percent FDR snp count cutoff must be true).')
#     parser.add_argument("-fdr", "--fdr", action="store", default=1e-6, type=float,\
#         help='SNP false discovery rate- based on simulation data with a 0.1 percent error rate (Q30)')
#     parser.add_argument("--min_fasta_reads", action="store", default=0, type=int,\
#         help='Minimum number of reads mapping to a scaffold to proceed with profiling it')
#     parser.add_argument("--scaffolds_to_profile", action="store",\
#         help='Path to a file containing a list of scaffolds to profile- if provided will ONLY profile those scaffolds')
#
#     # Read filtering cutoffs
#     parser.add_argument("-l", "--filter_cutoff", action="store", default=0.95, type=float, \
#         help='Minimum percent identity of read pairs to consensus to use the reads. Must be >, not >=')
#     parser.add_argument("--min_mapq", action="store", default=-1, type=int,\
#         help='Minimum mapq score of EITHER read in a pair to use that pair. Must be >, not >=')
#     parser.add_argument("--max_insert_relative", action="store", default=3, type=float, \
#         help='Multiplier to determine maximum insert size between two reads - default is to use 3x median insert size. Must be >, not >=')
#     parser.add_argument("--min_insert", action="store", default=50, type=int,\
#         help='Minimum insert size between two reads - default is 50 bp. If two reads are 50bp each and overlap completely, their insert will be 50. Must be >, not >=')
#
#     parser.add_argument('--store_everything', action='store_true', default=False,\
#         help="Store intermediate dictionaries in the pickle file; will result in significantly more RAM and disk usage")
#     parser.add_argument('--skip_mm_profiling', action='store_true', default=False,\
#         help="Dont perform analysis on an mm level; saves RAM and time")
#
#     # parser.add_argument("-g", "--genes", action="store", default=None, \
#     #     help='Optional genes file')
#     parser.add_argument('--debug', action='store_true', default=False, \
#         help ="Produce some extra output helpful for debugging")
#
#     # Parse
#     if (len(args) == 0 or args[0] == '-h' or args[0] == '--help'):
#         parser.print_help()
#         sys.exit(0)
#     else:
#         return parser.parse_args(args)
