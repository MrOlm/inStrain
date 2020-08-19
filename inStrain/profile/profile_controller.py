"""
This handles the argument parsing and logic for profiling bam files
"""

# Import packages
import os
import copy
import logging
from tqdm import tqdm
import multiprocessing

import inStrain.logUtils
import inStrain.SNVprofile
import inStrain.GeneProfile

import inStrain.profile.linkage
import inStrain.profile.profile_utilities
import inStrain.profile.snv_utilities


class BamProfileController(object):
    """
    Handle the logic of profiling
    """

    def __init__(self, bam, Fdb, sR2M, ISP_loc, **kwargs):
        """
        The requirements are bam, Fdb, and sR2M

        Arguments:
            bam = location of .bam file
            Fdb = dictionary listing fasta locations to profile
            sR2M = dictionary of scaffold -> read pair -> number of mm
            ISP_loc = location of ISP profile to store results
        """
        self.bam_loc = bam
        self.sR2M = sR2M
        self.Fdb = Fdb
        self.ISP_loc = ISP_loc
        self.kwargs = kwargs

    def main(self):
        """
        Profile the bam with inStrain using the options in kwargs
        """
        # Get arguments
        self.gen_prof_args()

        # Generate commands
        self.make_profile_commands()

        # Run the multiprocessing
        self.run_profile_processing()

        # Return the results as an SNVprofile
        return self.ISP

    def gen_prof_args(self):
        """
        Generate a set of arguments to be passed to prepare_comands
        """
        # Make a special copy of the arguments
        profArgs = copy.deepcopy(self.kwargs)
        scaff2sequence = profArgs.pop('s2s')
        self.scaff2sequence = scaff2sequence
        self.profArgs = profArgs
        self.scaffold_num = len(self.Fdb['scaffold'].unique())

        # Generate the null model
        fdr = profArgs.get('fdr', 1e-6)
        null_loc = os.path.dirname(__file__) + '/../helper_files/NullModel.txt'
        null_model = inStrain.profile.snv_utilities.generate_snp_model(null_loc, fdr=fdr)
        self.null_model = null_model

        # Parse for ProfileGenes as needed
        if self.kwargs.get('gene_file') is not None:
            self.gen_genes_args()
        else:
            self.profile_genes = False

    def gen_genes_args(self):
        """
        Parse arguments as needed to profile genes during the "merge" step
        """
        gene_file = self.kwargs.get('gene_file')

        inStrain.logUtils.log_checkpoint('Profile', 'Loading_genes', 'start')
        gene_database, gene2sequence = inStrain.GeneProfile.parse_genes(gene_file, **self.kwargs)
        self.profile_genes = True
        self.gene_database = gene_database
        self.gene2sequence = gene2sequence
        inStrain.logUtils.log_checkpoint('Profile', 'Loading_genes', 'end')

    def make_profile_commands(self):
        """
        Make a set of commands to run with profile
        """
        logging.debug('Creating commands')

        cmd_groups, Sprofile_dict, s2splits = prepare_commands(self.Fdb, self.bam_loc,
                                                               self.scaff2sequence, self.profArgs, self.sR2M)
        self.cmd_groups = cmd_groups
        self.Sprofile_dict = Sprofile_dict
        self.s2splits = s2splits

        logging.debug('There are {0} cmd groups'.format(len(cmd_groups)))

    def run_profile_processing(self):
        """
        Do the actual multiprocessing involved
        """
        # Establish command and result queues
        self.make_profile_queues()

        # Spawn worker processes for profiling splits
        self.spawn_profile_workers()

        # Get results from profiling splits
        self.recieve_profile_results()

        # Spwan worker processes for merging splits
        self.spawn_profile_merge_workers()

        # Get results from merging splits
        self.recieve_merge_results()

        # Collate and return results
        ISP = inStrain.profile.profile_utilities.gen_snv_profile(
            [s for s in self.Sprofiles if s is not None],
            ISP_loc=self.ISP_loc,
            **self.kwargs)

        # Store gene stuff if you have it
        if self.profile_genes:
            ISP.store('genes_fileloc', self.kwargs.get('gene_file'), 'value',
                      'Location of genes file that was used to call genes')
            ISP.store('genes_table', self.gene_database, 'pandas',
                      'Location of genes in the associated genes_file')
            if self.kwargs.get('store_everything', False):
                ISP.store('gene2sequence', self.gene2sequence, 'pickle',
                          'Dicitonary of gene -> nucleotide sequence')

        self.ISP = ISP

    def make_profile_queues(self):
        """
        Make the queues that will be used for profile multiprocessing
        """
        inStrain.logUtils.log_checkpoint("Profile", "initialize_multiprocessing", "start")

        # Get the arguments
        Sprofile_dict = self.Sprofile_dict
        cmd_groups = self.cmd_groups
        s2splits = self.s2splits

        # Create the queues
        ctx = multiprocessing.get_context('spawn')
        self.ctx = ctx

        manager = self.ctx.Manager()
        self.Sprofile_dict = manager.dict(Sprofile_dict)  # Holds a synced directory of splits

        self.log_list = ctx.Queue()
        self.split_cmd_queue = ctx.Queue()
        self.sprofile_cmd_queue = ctx.Queue()
        self.sprofile_result_queue = ctx.Queue()
        self.Sprofiles = []

        # Submit commands to the queues
        logging.debug("Submitting split tasks to queue")
        for i, cmd_group in enumerate(cmd_groups):
            self.split_cmd_queue.put(cmd_group)

        logging.debug("Submitting merge tasks to queue")
        for scaff, splits in s2splits.items():
            Sprofile = inStrain.profile.profile_utilities.ScaffoldSplitObject(splits)
            Sprofile.scaffold = scaff
            Sprofile.profile_genes = self.profile_genes

            if self.profile_genes:
                self.add_gene_info(Sprofile)

            self.sprofile_cmd_queue.put(Sprofile)

        inStrain.logUtils.log_checkpoint("Profile", "initialize_multiprocessing", "end")

    def add_gene_info(self, Sprofile):
        gdb = self.gene_database[self.gene_database['scaffold'] == Sprofile.scaffold]
        genes = set(gdb['gene'].tolist())
        if len(gdb) == 0:
            return
        else:
            Sprofile.gene_database = gdb
            Sprofile.gene2sequence = {g:s for g, s in self.gene2sequence.items() if g in genes}


    def spawn_profile_workers(self):
        """
        Spawn worker threads for profiling splits; or just run them all if a single core
        """
        p = int(self.kwargs.get('processes', 6))

        if p > 1:
            inStrain.logUtils.log_checkpoint("Profile", "SpawningSplitWorkers", "start")
            self.processes = []
            for i in range(0, p):
                self.processes.append(self.ctx.Process(
                    target=inStrain.profile.profile_utilities.split_profile_worker,
                    args=(self.split_cmd_queue,
                          self.Sprofile_dict,
                          self.log_list,
                          self.null_model,
                          self.bam_loc)))
            for proc in self.processes:
                proc.start()
            inStrain.logUtils.log_checkpoint("Profile", "SpawningSplitWorkers", "end")

        else:
            inStrain.profile.profile_utilities.split_profile_worker(self.split_cmd_queue,
                                                                    self.Sprofile_dict,
                                                                    self.log_list,
                                                                    self.null_model,
                                                                    self.bam_loc,
                                                                    single_thread=True)
            logging.info("Done profiling splits")

    def recieve_profile_results(self):
        """
        Get the results from the queues for profiling splits
        """
        p = int(self.kwargs.get('processes', 6))
        total_cmd_count = len(self.cmd_groups)

        if p > 1:
            # Set up progress bar
            pbar = tqdm(desc='Profiling splits: ', total=total_cmd_count)

            # Get the splits
            received_splits = 0
            while received_splits < total_cmd_count:
                try:
                    log = self.log_list.get()
                    logging.debug(log)
                    pbar.update(1)
                    received_splits += 1
                except KeyboardInterrupt:
                    for proc in self.processes:
                        proc.terminate()
                    break

            # Close progress bar
            pbar.close()
            inStrain.logUtils.log_checkpoint("Profile", "TerminatingSplitWorkers", "start")
            for proc in self.processes:
                proc.terminate()
            inStrain.logUtils.log_checkpoint("Profile", "TerminatingSplitWorkers", "end")

        else:
            # Get the splits
            received_splits = 0
            while received_splits < total_cmd_count:
                try:
                    log = self.log_list.get(timeout=5)
                    logging.debug(log)
                    received_splits += 1
                except:
                    logging.warning("Missing splits; {0} {1}".format(self.cmd_groups, self.split_cmd_queue))
                    assert False

    def spawn_profile_merge_workers(self):
        """
        Spawn worker threads for merging splits; or just run them all if a single core
        """
        p = int(self.kwargs.get('processes', 6))

        logging.debug('Establishing processes for merging')
        if p > 1:
            self.processes = []
            for i in range(0, p):
                self.processes.append(self.ctx.Process(
                    target=inStrain.profile.profile_utilities.merge_profile_worker,
                    args=(self.sprofile_cmd_queue,
                          self.Sprofile_dict,
                          self.sprofile_result_queue,
                          self.null_model)))
            for proc in self.processes:
                proc.start()

        else:
            inStrain.profile.profile_utilities.merge_profile_worker(
                self.sprofile_cmd_queue,
                self.Sprofile_dict,
                self.sprofile_result_queue,
                self.null_model,
                single_thread=True)

    def recieve_merge_results(self):
        """
        Get the results from queues after doing merges
        """
        p = int(self.kwargs.get('processes', 6))

        Sprofiles = []
        if p > 1:
            # Set up progress bar
            pbar = tqdm(desc='Merging splits and profiling genes: ', total=self.scaffold_num)

            # Get results
            received_profiles = 0
            while received_profiles < self.scaffold_num:
                try:
                    Sprofile = self.sprofile_result_queue.get()
                    if Sprofile is not None:
                        logging.debug(Sprofile.merge_log)
                        Sprofiles.append(Sprofile)
                    pbar.update(1)
                    received_profiles += 1
                except KeyboardInterrupt:
                    break

            # Close multi-processing
            for proc in self.processes:
                proc.terminate()

            # Close progress bar
            pbar.close()

        else:
            # Get results
            received_profiles = 0
            while received_profiles < self.scaffold_num:
                Sprofile = self.sprofile_result_queue.get(timeout=5)
                logging.debug(Sprofile.merge_log)
                Sprofiles.append(Sprofile)
                received_profiles += 1

        self.Sprofiles = Sprofiles


class profile_command():
    """
    This is a stupid object that just holds the arguments to profile a split
    """

    def __init__(self):
        pass


def prepare_commands(Fdb, bam, scaff2sequence, args, sR2M):
    """
    Make and iterate profiling commands
    Doing it in this way makes it use way less RAM
    """

    processes = args.get('processes', 6)
    s2p = args.get('s2p', None)
    if s2p is not None:
        SECONDS = min(60,
                      sum(calc_estimated_runtime(s2p[scaff]) for scaff in Fdb['scaffold'].unique()) / (processes + 1))
    else:
        SECONDS = 60

    cmd_groups = []
    Sdict = {}
    s2splits = {}

    seconds = 0
    cmds = []
    for scaff, db in Fdb.groupby('scaffold'):
        s2splits[scaff] = len(db)

        for i, row in db.iterrows():

            # make this command
            start = int(row['start'])
            end = int(row['end'])

            cmd = profile_command()
            cmd.scaffold = scaff
            cmd.R2M = sR2M[scaff]
            cmd.samfile = bam
            cmd.arguments = args
            cmd.start = start
            cmd.end = end
            cmd.split_number = int(row['split_number'])
            cmd.sequence = scaff2sequence[scaff][start:end + 1]

            # Add to the Sdict
            Sdict[scaff + '.' + str(row['split_number'])] = None

            # Add estimated seconds
            seconds += calc_estimated_runtime(s2p[scaff]) / s2splits[scaff]
            cmds.append(cmd)

            # See if you're done
            if seconds >= SECONDS:
                cmd_groups.append(cmds)
                seconds = 0
                cmds = []

    if len(cmds) > 0:
        cmd_groups.append(cmds)

    return cmd_groups, Sdict, s2splits


def calc_estimated_runtime(pairs):
    """
    Based on the number of mapped pairs, guess how long the split will take
    """
    SLOPE_CONSTANT = 0.0061401594694834305
    return (pairs * SLOPE_CONSTANT) + 0.2
