# Get the version
from ._version import __version__

# Import packages
import gc
import os
import sys
import h5py
import pysam
import logging
import argparse
import pandas as pd
from Bio import SeqIO
from subprocess import call
from datetime import datetime
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
import inStrain.logUtils

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

        self.shutdown(args)

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
        if args.run_statistics != None:
            inStrain.logUtils.process_logs(args.run_statistics)

    def shutdown(self, args):
        try:
            logloc = logging.getLoggerClass().root.handlers[0].baseFilename
        except:
            return
        logging.debug("inStrain complete; shutting down logger and printing run stats (log location = {0})".format(logloc))
        # reset_logging()
        logging.shutdown()
        # handlers = logging.getLogger('').handlers
        # for handler in handlers:
        #     handler.close()
        #     logging.getLogger('').removeHandler(handler)
        inStrain.logUtils.report_run_stats(logloc, most_recent=True, printToo=args.debug, debug=args.debug)

# def reset_logging():
#     '''
#     https://stackoverflow.com/questions/12034393/import-side-effects-on-logging-how-to-reset-the-logging-module
#     '''
#     manager = logging.root.manager
#     manager.disabled = logging.NOTSET
#     for logger in manager.loggerDict.values():
#         if isinstance(logger, logging.Logger):
#             logger.setLevel(logging.NOTSET)
#             logger.propagate = True
#             logger.disabled = False
#             logger.filters.clear()
#             handlers = logger.handlers.copy()
#             for handler in handlers:
#                 # Copied from `logging.shutdown`.
#                 try:
#                     handler.acquire()
#                     handler.flush()
#                     handler.close()
#                 except (OSError, ValueError):
#                     pass
#                 finally:
#                     handler.release()
#                 logger.removeHandler(handler)

class ProfileController():
    '''
    Main controller of the profile command
    '''

    def main(self, args):
        '''
        The main method when run on the command line
        '''
        # Parse arguments
        if args.output == 'inStrain':
            args.output = args.fasta.split(".")[0].split("/")[-1] #default prefix is now fasta prefix -alexcc 5/8/2019
        args = self.validate_arguments(args)
        message = """\
***************************************************
    ..:: inStrain profile Step 1. Filter reads ::..
***************************************************
        """
        logging.info(message)
        inStrain.logUtils.log_checkpoint("main_profile", "filter_reads", "start")

        # global s2l # make ths global so we can access it later.

        # Parse the args
        bam = args.bam
        vargs = vars(args)
        del vargs['bam']

        # Set up .fasta file
        inStrain.logUtils.log_checkpoint("FilterReads", "load_fasta", "start")
        FAdb, s2s = self.load_fasta(args)
        s2l = {s:len(s2s[s]) for s in list(s2s.keys())}
        inStrain.logUtils.log_checkpoint("FilterReads", "load_fasta", "end")

        # Load dictionary of paired reads
        scaffolds = list(FAdb['scaffold'].unique())

        detailed_report = vargs.get('detailed_mapping_info', False)

        if detailed_report:
            Rdic, RR, dRR = inStrain.filter_reads.load_paired_reads(bam, scaffolds, **vargs)
        else:
            Rdic, RR = inStrain.filter_reads.load_paired_reads(bam, scaffolds, **vargs)
            dRR = None

        logging.info("{0:,} read pairs remain after filtering".format(RR['filtered_pairs'].tolist()[0]))

        if RR['filtered_pairs'].tolist()[0] == 0:
            logging.error("Because no read pairs remain I'm going to crash now. Maybe this is failing because you dont have paired reads (in which case you should adjust --pairing_filter option), or maybe its failing because the mapper you used uses full fasta headers (in which case you should use the flag --use_full_fasta_header)")
            raise Exception('No paired reads detected; see above message and log')

        # Get scaffold to paired reads (useful for multiprocessing)
        s2p = RR.set_index('scaffold')['filtered_pairs'].to_dict()

        # Get the read length
        rl = float(RR.loc[0, 'mean_pair_length'])

        # Filter the .fasta file
        FAdb = self.filter_fasta(FAdb, s2p, s2l, rl, **vargs)
        assert len(FAdb) > 0, "No scaffolds passed initial filtering based on numbers of mapped reads"

        if args.skip_mm_profiling:
            newRdic = {}
            for s, p2i in Rdic.items():
                newRdic[s] = set(p2i.keys())
            del Rdic
            Rdic = newRdic

        # if args.debug:
        #     for att in ['Rdic', 'RR', 'FAdb', 's2s']:
        #         logging.debug("RAM PROFILE: reads {0} {1:.2f} Mb".format(att,
        #             sys.getsizeof(eval(att))/1e6))

        # Profile the .bam file
        vargs['s2s'] = s2s
        vargs['s2p'] = s2p

        inStrain.logUtils.log_checkpoint("main_profile", "filter_reads", "end")
        message = """\
***************************************************
.:: inStrain profile Step 2. Profile scaffolds ::..
***************************************************
        """
        logging.info(message)
        inStrain.logUtils.log_checkpoint("main_profile", "profile_scaffolds", "start")

        Sprofile = inStrain.profileUtilities.profile_bam(bam, FAdb, Rdic, **vargs)

        # Add the read report
        Sprofile.store('mapping_info', RR, 'pandas', "Report on reads")
        if dRR != None:
            Sprofile.store('detailed_mapping_info', dRR, 'pandas', "Details report on reads")
            del dRR

        logging.debug("Storing Rdic")
        if args.skip_mm_profiling:
            Sprofile.store('Rdic', Rdic, 'pickle', 'list of filtered read pairs')
        else:
            Sprofile.store('Rdic', Rdic, 'dictionary', 'Read pair -> mismatches')
        logging.debug("Done storing Rdic")

        # Store the .fasta location
        Sprofile.store('fasta_loc', os.path.abspath(args.fasta), 'value', 'Location of .fasta file used during profile')
        Sprofile.store('scaffold2length', s2l, 'dictionary', 'Dictionary of scaffold 2 length')

        # Save
        logging.info('Storing output')
        self.write_output(Sprofile, args)

        # Run the rest of things
        args.IS = Sprofile.location

        inStrain.logUtils.log_checkpoint("main_profile", "profile_scaffolds", "end")
        # See if you can profile genes as well
        message = """\
***************************************************
  .:: inStrain profile Step 3. Profile genes ::..
***************************************************
        """
        logging.info(message)
        if args.gene_file != None:
            inStrain.logUtils.log_checkpoint("main_profile", "profile_genes", "start")
            args.IS = Sprofile.location
            Controller().profile_genes_operation(args)
            inStrain.logUtils.log_checkpoint("main_profile", "profile_genes", "end")
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
            inStrain.logUtils.log_checkpoint("main_profile", "genome_wide", "start")
            args.IS = Sprofile.location
            Controller().genome_wide_operation(args)
            inStrain.logUtils.log_checkpoint("main_profile", "genome_wide", "end")
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
            inStrain.logUtils.log_checkpoint("main_profile", "making_plots", "start")
            args.IS = Sprofile.location
            args.plots = 'a'
            Controller().plot_operation(args)
            inStrain.logUtils.log_checkpoint("main_profile", "making_plots", "end")
        else:
            logging.info('Nevermind! You chose to skip making plots')

        # Final message
        message = """\
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ..:: inStrain profile finished ::..

Output tables........ {0}
Figures.............. {1}
Logging.............. {2}

See documentation for output descriptions - https://instrain.readthedocs.io/en/latest/

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        """.format(Sprofile.get_location('output'), \
            Sprofile.get_location('figures'),
            Sprofile.get_location('log'))
        logging.info(message)

        return Sprofile

    def load_fasta(self, args):
        '''
        Load the sequences to be profiled

        Return a table listing scaffold name, start, end
        '''
        if args.use_full_fasta_header:
            #scaff2sequence = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"), key_function=_get_description)
            scaff2sequence = {r.description:r.seq.upper() for r in SeqIO.parse(args.fasta, "fasta")}
        else:
            #scaff2sequence = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
            scaff2sequence = {r.id:r.seq.upper() for r in SeqIO.parse(args.fasta, "fasta")}

        # Get scaffold2length
        s2l = {s:len(scaff2sequence[s]) for s in list(scaff2sequence.keys())}

        # Generate splits
        table = defaultdict(list)
        WINDOW_LEN = args.window_length
        for scaffold, sLen in s2l.items():
            for i, (split_start, split_end) in enumerate(iterate_splits(sLen, WINDOW_LEN)):
                table['scaffold'].append(scaffold)
                table['split_number'].append(i)
                table['start'].append(split_start)
                table['end'].append(split_end)
        Fdb = pd.DataFrame(table)

        # This just takes too long, sadly
        # _validate_splits(Fdb, s2l)

        if args.scaffolds_to_profile != None:
            Fdb = Fdb[Fdb['scaffold'].isin(args.scaffolds_to_profile)]
            s2s = {scaff:seq for scaff, seq in scaff2sequence.items() if scaff in args.scaffolds_to_profile}

        if len(Fdb) == 0:
            logging.error("The provided scaffold list has no overlap with the provided .fasta file!")
            logging.error("Example scaffolds in list: {0}".format("\n".join(args.scaffolds_to_profile)))
            sys.exit()

        return Fdb, scaff2sequence # also return s2l - alexcc 5/9/2019: Nah, make it scaff2sequence (s2s) (M.O. 6/10/19)

    # def load_paired_reads(self, args, FAdb):
    #     '''
    #     Load paired reads to be profiled
    #
    #     Return a dictionary of read pair -> info
    #     '''
    #     # Load paired reads
    #     scaffolds = list(FAdb['scaffold'].unique())
    #     scaff2pair2info, scaff2total = inStrain.filter_reads.get_paired_reads_multi(args.bam, scaffolds, processes=args.processes, ret_total=True)
    #
    #     # Merge into single
    #     pair2info = {}
    #     for scaff, p2i in scaff2pair2info.items():
    #         for p, i in p2i.items():
    #             pair2info[p] = i
    #
    #     # Make read report
    #     logging.info('Making read report')
    #     RR = inStrain.filter_reads.makeFilterReport(scaff2pair2info, scaff2total, pair2info=pair2info, **vars(args))
    #
    #     # Filter the dictionary
    #     logging.info('Filtering reads')
    #     pair2infoF = inStrain.filter_reads.filter_paired_reads_dict(pair2info,
    #         **vars(args))
    #
    #     # Make a report on these pair reads
    #     #make_mapping_info(pair2info, pair2infoF, args)
    #
    #     return pair2infoF, RR

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
        args.bam = prepare_bam_fie(args)

        # Load the list of scaffolds
        args.scaffolds_to_profile = load_scaff_list(args.scaffolds_to_profile)

        # Fix the fdr
        if args.fdr == 0:
            args.fdr = 1e-6

        # Handle database mode
        if args.database_mode:
            args.min_read_ani = 0.92
            args.skip_mm_profiling = True
            args.min_genome_coverage = 1

        # Make sure you have a .stb if needed
        if args.min_genome_coverage != 0:
            assert args.stb != [], 'If you adjust the minimum genome coverage, you need to provide an .stb!'

        return args

    def write_output(self, Sprofile, args):
        '''
        Write output files
        '''
        logging.debug("Writing output files now")

        for t in ['SNVs', 'scaffold_info', 'SNVs', 'linkage']:
            Sprofile.generate(t)
        Sprofile.generate('mapping_info', **vars(args))

        # out_base = Sprofile.get_location('output') + os.path.basename(Sprofile.get('location')) + '_'
        #
        # # Write the scaffold profile
        # db = Sprofile.get_nonredundant_scaffold_table()
        # db.to_csv(out_base + 'scaffold_info.tsv', index=False, sep='\t')
        #
        # # Write the SNP frequencies
        # db = Sprofile.get_nonredundant_snv_table()
        # db.to_csv(out_base + 'SNVs.tsv', index=False, sep='\t')
        #
        # # Write the linkage
        # db = Sprofile.get_nonredundant_linkage_table()
        # db.to_csv(out_base + 'linkage.tsv', index=False, sep='\t')
        #
        # # Write the read report
        # RR = Sprofile.get('mapping_info')
        # inStrain.filter_reads.write_mapping_info(RR, out_base + 'mapping_info.tsv', **vars(args))

    def filter_fasta(self, FAdb, s2p, s2l, rl, **kwargs):
        '''
        Filter the .fasta file based on the min number of mapped paired reads
        '''
        min_reads = kwargs.get('min_scaffold_reads', 0)
        min_genome_coverage = kwargs.get('min_genome_coverage', 0)

        if min_genome_coverage > 0:
            FAdb = _filter_genome_coverage(FAdb, s2l, s2p, rl, min_genome_coverage, kwargs.get('stb'))

        if len(FAdb) == 0:
            return FAdb

        # Remove scaffolds without the min number of reads
        FAdb = FAdb[[True if (s2p[s] >= min_reads) else False for s in FAdb['scaffold']]]

        # Sort scaffolds based on the number of reads
        FAdb['filtered_pairs'] = FAdb['scaffold'].map(s2p)
        FAdb = FAdb.sort_values('filtered_pairs', ascending=False)

        return FAdb

def _filter_genome_coverage(FAdb, s2l, s2p, rl, min_genome_coverage, stb_loc):
    '''
    Calcualte the coverage of genomes based on the read filtering, and only keep scaffolds that are above the threshold

    stb_loc should be a list, direct from argument parser
    '''
    cdb = FAdb.drop_duplicates(subset=['scaffold'])
    cdb['read_pairs'] = cdb['scaffold'].map(s2p)
    cdb['length'] = cdb['scaffold'].map(s2l)

    stb = inStrain.genomeUtilities.load_scaff2bin(stb_loc)
    cdb = inStrain.genomeUtilities._add_stb(cdb, stb)

    xdb = cdb.groupby('genome')[['read_pairs', 'length']].agg(sum).reset_index()
    xdb['genome_bases'] = xdb['read_pairs'] * rl
    xdb['coverage'] = xdb['genome_bases'] / xdb['length']
    genome_to_rm = set(xdb[xdb['coverage'] < min_genome_coverage]['genome'].tolist())

    scaffolds_to_rm_1 = set(cdb[cdb['genome'].isin(genome_to_rm)]['scaffold'].tolist())
    scaffolds_to_rm_2 = set(cdb[cdb['genome'].isna()]['scaffold'].tolist())
    scaffolds_to_rm = scaffolds_to_rm_1.union(scaffolds_to_rm_2)

    logging.info("{0} of {1} genomes have less than {2}x estimated coverage".format(
            len(genome_to_rm), len(xdb['genome'].unique()), min_genome_coverage))
    logging.info("{0} of the original {1} scaffolds are removed ({2} have a low coverage genome; {3} have no genome)".format(
            len(scaffolds_to_rm), len(cdb['scaffold'].unique()), len(scaffolds_to_rm_1), len(scaffolds_to_rm_2)))

    return FAdb[~FAdb['scaffold'].isin(scaffolds_to_rm)]

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

def iterate_splits(sLen, WINDOW_LEN):
    '''
    Splits are 0-based and double-inclusive
    '''
    numberChunks = sLen // WINDOW_LEN + 1
    chunkLen = int(sLen / numberChunks)

    #print("Scaffold length of {0}, window length of {1}, {2} splits of {3}".format(sLen, WINDOW_LEN, numberChunks, chunkLen))

    start = 0
    end = 0
    for i in range(numberChunks):
        if i + 1 == numberChunks:
            yield start, sLen - 1
        else:
            end += chunkLen
            yield start, end - 1
            start += chunkLen

def _validate_splits(Fdb, s2l):
    '''
    Splits are 0-based and double-inclusive
    '''
    for scaffold, db in Fdb.groupby('scaffold'):
        db['len'] = db['end'] - db['start'] + 1
        if db['len'].sum() != s2l[scaffold]:
            print(db)
        assert db['len'].sum() == s2l[scaffold], [db['len'].sum(), s2l[scaffold]]
        assert db['start'].min() == 0
        assert db['end'].max() == s2l[scaffold] - 1


def setup_logger(loc):
    ''' set up logger such that DEBUG goes only to file, rest go to file and console '''

    # Cancel if a logger already exists:
    if logging.getLogger('').handlers:
        return

    # set up logging everything to file
    logging.basicConfig(level=logging.DEBUG,
                       format='%(asctime)s %(levelname)-8s %(message)s',
                       datefmt='%y-%m-%d %H:%M:%S',
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

def report_run_stats(logloc, most_recent=True, printToo=True):
    if logloc == None:
        return

    Ldb = load_log(logloc)
    print(Ldb)

def load_log(logfile):
    table = defaultdict(list)
    with open(logfile) as o:
        prev_line = None
        for line in o.readlines():
            line = line.strip()

            # load new inStrain run
            if 'inStrain version' in line:
                linewords = [x.strip() for x in line.split()]
                epoch_time = log_fmt_to_epoch("{0} {1}".format(linewords[0], linewords[1]))

                table['log_type'].append('program_start')
                table['time'].append(epoch_time)
                table['parsable_string'].append("version={0}".format(linewords[5]))

            # Load  profile RAM and multiprocessing reporting
            elif 'RAM. System has' in line:
                linewords = [x.strip() for x in line.split()]
                pstring = "scaffold={0};PID={1};status={2};process_RAM={3};system_RAM={4};total_RAM={5}".format(
                            linewords[0], linewords[2], linewords[3], linewords[7], linewords[11], linewords[13])

                table['log_type'].append('Profile_PID_RAM')
                table['time'].append(linewords[5])
                table['parsable_string'].append(pstring)
                # table['scaffold'].append(linewords[0])
                # table['PID'].append(linewords[2])
                # table['status'].append(linewords[3])
                # table['time'].append(linewords[5])
                # table['process_RAM'].append(linewords[7])
                # table['system_RAM'].append(linewords[11])
                # table['total_RAM'].append(linewords[13])
            prev_line = line

    Ldb = pd.DataFrame(table)
    return Ldb

def log_fmt_to_epoch(ttime):
    oldformat = '%m-%d %H:%M'
    print(ttime)
    datetimeobject = datetime.strptime(ttime,oldformat)
    print(datetimeobject)
    return datetimeobject.timestamp()

def prepare_bam_fie(args):
    '''
    Make this a .bam file
    '''
    bam = args.bam

    if bam[-4:] == '.sam':
        logging.info("You gave me a sam- I'm going to make it a .bam now")

        bam = _sam_to_bam(bam)
        bam = _sort_index_bam(bam)

    elif bam[-4:] == '.bam':
        # If there's an index, assume its sorted
        if (os.path.exists(bam + '.bai')) | ((os.path.exists(bam[:-4] + '.bai'))):
            pass

        # If there's not an index...
        else:
            bam = _sort_index_bam(bam, rm_ori=False)

    if os.stat(bam).st_size == 0:
        logging.error("Failed to generated a sorted .bam file! Make sure you have "+\
            "samtools version 1.6 or greater.")
        sys.exit()

    # Do a quick sanity check
    try:
         pysam.AlignmentFile(bam).mapped
    except ValueError:
        logging.error("It seems like the .bam file could not be indexed!"  +\
                    "Make sure you have samtools version 1.6 or greater.")
        sys.exit()

    return bam

def _sam_to_bam(sam):
    '''
    From the location of a .sam file, convert it to a bam file and retun the location
    '''
    if sam[-4:] != '.sam':
        print('Sam file needs to end in .sam')
        sys.exit()

    bam = sam[:-4] + '.bam'
    logging.info("Converting {0} to {1}".format(sam, bam))
    cmd = ['samtools', 'view','-S','-b', sam, '>', bam]
    print(' '.join(cmd))
    call(' '.join(cmd), shell=True)

    return bam

def _sort_index_bam(bam, rm_ori=False):
    '''
    From a .bam file, sort and index it. Remove original if rm_ori
    Return path of sorted and indexed bam
    '''
    if bam[-4:] != '.bam':
        logging.error('Bam file needs to end in .bam')
        sys.exit()

    if 'sorted.bam' not in bam:
        logging.info("sorting {0}".format(bam))
        sorted_bam = bam[:-4] + '.sorted.bam'
        cmd = ['samtools', 'sort', bam, '-o', sorted_bam]
        print(' '.join(cmd))
        call(cmd)
    else:
        sorted_bam = bam
        rm_ori = False

    logging.info("Indexing {0}".format(sorted_bam))
    cmd = ['samtools', 'index', sorted_bam, sorted_bam + '.bai']
    print(' '.join(cmd))
    call(cmd)

    if rm_ori:
        logging.info("Deleting {0}".format(bam))
        os.remove(bam)

    return sorted_bam
    #return convert_table(Ldb)
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
#     parser.add_argument("--min_scaffold_reads", action="store", default=0, type=int,\
#         help='Minimum number of reads mapping to a scaffold to proceed with profiling it')
#     parser.add_argument("--scaffolds_to_profile", action="store",\
#         help='Path to a file containing a list of scaffolds to profile- if provided will ONLY profile those scaffolds')
#
#     # Read filtering cutoffs
#     parser.add_argument("-l", "--min_read_ani", action="store", default=0.95, type=float, \
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
