# Get the version
from ._version import __version__

# Import packages
import gc
import os
import sys
import h5py
import copy
import pysam
import logging
import argparse
import pandas as pd
pd.set_option('mode.chained_assignment',None)
import functools
import logging
import struct
import sys

from subprocess import call
from datetime import datetime
from collections import defaultdict

# Import inStrain stuff
import inStrain.profile
import inStrain.profile.samtools_ops

import inStrain.plotting
import inStrain.plotting.plotting_controller

import inStrain.filter_reads
import inStrain.utils
import inStrain.readComparer
import inStrain.GeneProfile
import inStrain.genomeUtilities
import inStrain.quickProfile
import inStrain.SNVprofile
import inStrain.logUtils
import inStrain.compare_controller

class Controller():
    '''
    Controller of the whole shebang
    '''

    def main(self, args):
        ''' Parse user options and call the correct pipeline'''
        # Patch multiprocessing
        patch_mp_connection_bpo_17560()

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

        if args.operation == "check_deps":
            self.checkdep_operation(args)

        self.shutdown(args)

    def profile_operation(self, args):
        ProfileController(args).main()

    def compare_operation(self, args):
        inStrain.compare_controller.CompareController(args).main()

    def filter_reads_operation(self, args):
        inStrain.filter_reads.Controller().main(args)

    def profile_genes_operation(self, args):
        inStrain.GeneProfile.Controller().main(args)

    def genome_wide_operation(self, args):
        inStrain.genomeUtilities.Controller().main(args)

    def quick_profile_operation(self, args):
        inStrain.quickProfile.main(args)

    def plot_operation(self, args):
        inStrain.plotting.plotting_controller.PlottingController(args).main()
        #inStrain.plottingUtilities.main(args)

    def other_operation(self, args):
        # Check if you should convert IS profile
        if args.old_IS != None:
            inStrain.SNVprofile.convert_SNVprofile(args.old_IS)
        if args.run_statistics != None:
            inStrain.logUtils.process_logs(args.run_statistics)

    def checkdep_operation(self, args):
        message = inStrain.utils.gen_dependency_report()
        print(message)

    def shutdown(self, args):
        try:
            logloc = logging.getLoggerClass().root.handlers[0].baseFilename
        except:
            return
        logging.debug("inStrain complete; shutting down logger and printing run stats (log location = {0})".format(logloc))
        logging.shutdown()
        inStrain.logUtils.report_run_stats(logloc, most_recent=True, printToo=args.debug, debug=args.debug)

class ProfileController(object):
    '''
    Main controller of the profile command
    '''
    def __init__(self, args):
        '''
        Set all of the command line arguments in the "args" attribute

        Doing it this way lets your pass the arguments to other controllers
        '''
        self.args = args
        self.ori_args = copy.deepcopy(args)

    def main(self):
        '''
        The main method when run on the command line
        '''
        # Parse arguments
        self.validate_arguments()

        # Filter reads
        self.profile_filter_reads()

        # Profile
        self.run_profile()

        # Profile genes
        # self.profile_profile_genes()

        # Make things genome_wide
        self.profile_genome_wide()

        # Make plots
        self.profile_plots()

        # Final message
        self.write_final_message()

        return self.ISP

    def validate_arguments(self):
        '''
        Do some parsing, start up a logger
        '''
        # Get out the "args" to manipulate it
        args = self.args

        # default prefix is now fasta prefix -alexcc 5/8/2019
        if args.output == 'inStrain':
            args.output = args.fasta.split(".")[0].split("/")[-1]

        # Set up "base"
        out_base = args.output

        # Set up Logger
        outbase = out_base
        ISP = inStrain.SNVprofile.SNVprofile(outbase)
        log_loc = ISP.get_location('log') + 'log.log'
        setup_logger(log_loc)

        # Make the bam file if you need to; remove it from args
        self.bam = inStrain.profile.samtools_ops.prepare_bam_fie(args.bam, args.processes)
        del self.args.bam

        # Load the list of scaffolds
        args.scaffolds_to_profile = inStrain.profile.fasta.load_scaff_list(
                                                    args.scaffolds_to_profile)

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

        # See if you have genes
        if args.gene_file != None:
            self.genes_file_present = True

        self.args = args
        self.args.IS = ISP.get('location')
        self.ISP = ISP
        self.kwargs = vars(self.args)

    def profile_filter_reads(self):
        '''
        Call the filter reads module as run with "profile"
        '''
        message = """\
***************************************************
    ..:: inStrain profile Step 1. Filter reads ::..
***************************************************
        """
        logging.info(message)
        inStrain.logUtils.log_checkpoint("main_profile", "filter_reads", "start")

        # Profile reads
        Rdic, RR, fasta_db, scaff2sequence, s2l = \
                inStrain.filter_reads.Controller().main_from_profile(self.ISP,
                                                    self.bam, **self.kwargs)

        # Store results
        self.RR = RR
        self.Rdic = Rdic
        self.fasta_db = fasta_db
        self.scaff2sequence = scaff2sequence
        self.scaffold2length = s2l

        # Parse results
        inStrain.logUtils.log_checkpoint("FilterReads", "parse_results", "start")
        self.parse_filter_reads()
        inStrain.logUtils.log_checkpoint("FilterReads", "parse_results", "end")

        inStrain.logUtils.log_checkpoint("main_profile", "filter_reads", "end")

    def parse_filter_reads(self):
        '''
        Parse filter read results in the context of profile
        '''
        # Parse results
        self.scaffold2pairs = self.RR.set_index('scaffold')['filtered_pairs'].to_dict()
        self.readlength = float(self.RR.loc[0, 'mean_pair_length'])

        # Filter the .fasta file with these results
        self.fasta_db = inStrain.profile.fasta.filter_fasta(self.fasta_db,
                        self.scaffold2pairs, self.scaffold2length,
                        self.readlength, **self.kwargs)

        # Change up Rdic if needed
        if self.args.skip_mm_profiling:
            newRdic = {}
            for s, p2i in self.Rdic.items():
                newRdic[s] = set(p2i.keys())
            self.Rdic = newRdic
            self.ISP.store('Rdic', self.Rdic, 'pickle', 'list of filtered read pairs')
        else:
            self.ISP.store('Rdic', self.Rdic, 'dictionary', 'Read pair -> mismatches')


        # Store results
        self.ISP.store('mapping_info', self.RR, 'pandas', "Report on reads")
        self.ISP.store('fasta_loc', os.path.abspath(self.args.fasta), 'value', 'Location of .fasta file used during profile')
        self.ISP.store('scaffold2length', self.scaffold2length, 'dictionary', 'Dictionary of scaffold 2 length')

        # Print status report
        unfiltered_pairs = self.RR['unfiltered_pairs'].iloc[0]
        filterd_pairs = self.RR['filtered_pairs'].iloc[0]
        pair_length = self.readlength

        status = ''
        status += "{0:.1f}% of reads were removed during filtering\n".format(
                    ((unfiltered_pairs - filterd_pairs) / unfiltered_pairs)*100)
        status += "{0:,} read pairs remain ({1:#.4g} Gbp)".format(
                    filterd_pairs, (filterd_pairs * pair_length)/1e9)
        logging.info(status)

        # Handle exceptions
        if self.RR['filtered_pairs'].tolist()[0] == 0:
            logging.error(
                "Because no read pairs remain I'm going to crash now. Maybe this is failing because you dont have paired reads (in which case you should adjust --pairing_filter option), or maybe its failing because the mapper you used uses full fasta headers (in which case you should use the flag --use_full_fasta_header)")
            raise Exception('No paired reads detected; see above message and log')

        if len(self.fasta_db) <= 0:
            logging.error("No scaffolds passed initial filtering based on numbers of mapped reads")
            raise Exception('No scaffolds detected; see above message and log')

    def run_profile(self):
        '''
        Call the actual profile module
        '''
        message = """\
***************************************************
.:: inStrain profile Step 2. Profile scaffolds ::..
***************************************************
        """
        logging.info(message)
        inStrain.logUtils.log_checkpoint("main_profile", "profile_scaffolds", "start")

        # Do some argument handling
        self.kwargs['s2s'] = self.scaff2sequence
        self.kwargs['s2p'] = self.scaffold2pairs

        # Call the module
        self.ISP = inStrain.profile.profile_bam(self.bam, self.fasta_db,
                                                self.Rdic, self.ISP.get('location'),
                                                **self.kwargs)

        # Write output
        inStrain.logUtils.log_checkpoint("Profile", "store_output", "start")
        self.write_output()
        inStrain.logUtils.log_checkpoint("Profile", "store_output", "end")

        inStrain.logUtils.log_checkpoint("main_profile", "profile_scaffolds", "end")

    def write_output(self):
        '''
        Write output files
        '''
        logging.debug("Writing output files now")

        for t in ['SNVs', 'scaffold_info', 'SNVs', 'linkage', 'gene_info']:
            self.ISP.generate(t)
        self.ISP.generate('mapping_info', **self.kwargs)

    def profile_profile_genes(self):
        '''
        Call profile genes from the "profile" module
        '''
        message = """\
***************************************************
  .:: inStrain profile Step 3. Profile genes ::..
***************************************************
        """
        logging.info(message)
        args = self.args

        if args.gene_file != None:
            inStrain.logUtils.log_checkpoint("main_profile", "profile_genes", "start")
            Controller().profile_genes_operation(copy.deepcopy(args))
            inStrain.logUtils.log_checkpoint("main_profile", "profile_genes", "end")
        else:
            logging.info('Nevermind! You didnt include a genes file')

    def profile_genome_wide(self):
        '''
        Call genome_wide from "profile" module
        '''
        message = """\
***************************************************
.:: inStrain profile Step 4. Make genome-wide ::..
***************************************************
        """
        logging.info(message)
        args = self.args

        if not args.skip_genome_wide:
            inStrain.logUtils.log_checkpoint("main_profile", "genome_wide", "start")
            Controller().genome_wide_operation(copy.deepcopy(args))
            inStrain.logUtils.log_checkpoint("main_profile", "genome_wide", "end")
        else:
            logging.info('Nevermind! You chose to skip genome_wide')

    def profile_plots(self):
        '''
        Call plotting function from "profile" module
        '''
        # Generate plots
        message = """\
***************************************************
 .:: inStrain profile Step 5. Generate plots ::..
***************************************************
        """
        logging.info(message)
        args = self.args

        if not args.skip_plot_generation:
            inStrain.logUtils.log_checkpoint("main_profile", "making_plots", "start")
            args.plots = 'a'
            Controller().plot_operation(args)
            inStrain.logUtils.log_checkpoint("main_profile", "making_plots", "end")
        else:
            logging.info('Nevermind! You chose to skip making plots')

    def write_final_message(self):
        Sprofile = self.ISP
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

    try:
        logging.debug(inStrain.utils.gen_dependency_report())
    except:
        logging.warning('Failed to log dependencies')


def patch_mp_connection_bpo_17560():
    """
    From https://stackoverflow.com/questions/47776486/python-struct-error-i-format-requires-2147483648-number-2147483647

    Apply PR-10305 / bpo-17560 connection send/receive max size update

    See the original issue at https://bugs.python.org/issue17560 and
    https://github.com/python/cpython/pull/10305 for the pull request.

    This only supports Python versions 3.3 - 3.7, this function
    does nothing for Python versions outside of that range.

    """
    patchname = "Multiprocessing connection patch for bpo-17560"
    if not (3, 3) < sys.version_info < (3, 8):
        # print(
        #     patchname + " not applied, not an applicable Python version: %s",
        #     sys.version
        # )
        return

    from multiprocessing.connection import Connection

    orig_send_bytes = Connection._send_bytes
    orig_recv_bytes = Connection._recv_bytes
    if (
        orig_send_bytes.__code__.co_filename == __file__
        and orig_recv_bytes.__code__.co_filename == __file__
    ):
        #print(patchname + " already applied, skipping")
        return

    @functools.wraps(orig_send_bytes)
    def send_bytes(self, buf):
        n = len(buf)
        if n > 0x7fffffff:
            pre_header = struct.pack("!i", -1)
            header = struct.pack("!Q", n)
            self._send(pre_header)
            self._send(header)
            self._send(buf)
        else:
            orig_send_bytes(self, buf)

    @functools.wraps(orig_recv_bytes)
    def recv_bytes(self, maxsize=None):
        buf = self._recv(4)
        size, = struct.unpack("!i", buf.getvalue())
        if size == -1:
            buf = self._recv(8)
            size, = struct.unpack("!Q", buf.getvalue())
        if maxsize is not None and size > maxsize:
            return None
        return self._recv(size)

    Connection._send_bytes = send_bytes
    Connection._recv_bytes = recv_bytes

   # print(patchname + " applied")
