"""
This handles the argument parsing and logic for profiling bam files
"""

# Import packages
import os
import copy
import logging
import pandas as pd
from tqdm import tqdm
import multiprocessing
from collections import defaultdict

import inStrain.readComparer
import inStrain.compare_utils

class CompareController(object):
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
        self.kwargs = vars(self.args)

        # Make null model for SNP calling
        fdr = self.kwargs.get('fdr', 1e-6)
        null_loc = os.path.dirname(__file__) + '/helper_files/NullModel.txt'
        self.null_model = inStrain.profile.snv_utilities.generate_snp_model(null_loc, fdr=fdr)

    def main(self):
        '''
        The main method when run on the command line
        '''
        self.parse_arguments()

        # Set up which scaffolds will be compared when
        self.create_scaffoldcomparison_objects()

        # Actually do the processing
        self.run_comparisons()

        # Store the results
        self.store_results()

    def parse_arguments(self):
        """
        Parse the arguments and add them this object as attributes
        """
        args = self.args

        # Set up output object and log
        outbase = args.output
        RCprof = inStrain.SNVprofile.SNVprofile(outbase)
        log_loc = RCprof.get_location('log') + 'log.log'
        inStrain.controller.setup_logger(log_loc)
        self.RCprof = RCprof

        # Set up list of input IS profiles
        inputs = list(args.input)
        assert len(inputs) > 1, "You need to have more than one input .IS file"
        self.inputs = inputs

        # Set up a list of scaffolds to compare
        if args.scaffolds is not None:
            self.scaffolds_to_compare = inStrain.profile.fasta.load_scaff_list(args.scaffolds)
        else:
            self.scaffolds_to_compare = None

        # Load the scaffold to bin file
        if args.stb is not None:
            if len(args.stb) > 0:
                # Load it and integrate with above
                print(args.stb)
                assert False, "MATT, YOU NEED TO SUPPORT STB FILES!"

    def create_scaffoldcomparison_objects(self):
        """
        Figure out which scaffolds need to be compared
        """
        # Load scaffolds to compare
        if self.scaffolds_to_compare is not None:
            scaffolds_to_compare = set(self.scaffolds_to_compare)
        else:
            scaffolds_to_compare = None

        # Load ScaffoldComparison objects
        valid_SCs, scaffold2length = make_scaffoldcomparison_objects(self.inputs, scaffolds_to_compare)
        self.SC_objects = valid_SCs
        self.scaffold2length = scaffold2length

        # Establish ScaffoldComparison groups
        group_length = 10000000

        #Cdb = calc_scaff_sim_matrix(valid_SCs)
        #SC_groups = establish_SC_groups(valid_SCs, Cdb, group_length)
        SC_groups = group_Scaffold_objects(valid_SCs, group_length)

        self.scaffold_comparison_groups = SC_groups

    def run_comparisons(self):
        """
        Run each group one at a time, storing the results in the meantime
        """
        cdbs = [] # Store scaffold comparison information
        mdbs = [] # Store mismatch locations
        pair2mm2covOverlaps = [] # Store coverage overlap locations
        order = [] # Store order of scaffolds

        groups = len(self.scaffold_comparison_groups)
        for i, SCgroup in enumerate(self.scaffold_comparison_groups):
            logging.info(f'Running group {i} of {groups}')
            SCgroup.load_cache()
            results = inStrain.compare_utils.run_compare_multiprocessing(SCgroup.cmd_queue, SCgroup.result_queue,
                                                                         self.null_model, num_to_run=len(SCgroup.scaffolds),
                                                                         **self.kwargs)
            for result in results:
                if result is not None:
                    Cdb, Mdb, pair2mm2covOverlap, scaffold = result
                    for item, lis in zip([Cdb, Mdb, pair2mm2covOverlap, scaffold], [cdbs, mdbs, pair2mm2covOverlaps, order]):
                        lis.append(item)

            SCgroup.purge_cache()

        # Process results
        self.process_results(cdbs, mdbs, pair2mm2covOverlaps, order)

    def process_results(self, cdbs, mdbs, pair2mm2covOverlaps, order):
        """
        Merge and store results
        """
        if len(pair2mm2covOverlaps) > 0:
            scaff2pair2mm2overlap = {}
            for scaff, pair2mm2covOverlap in zip(order, pair2mm2covOverlaps):
                scaff2pair2mm2overlap[scaff] = pair2mm2covOverlap
        else:
            scaff2pair2mm2overlap = None

        if len(mdbs) > 0:
            Mdb = pd.concat(mdbs, sort=False)
        else:
            Mdb = None

        # Do some typing of this DataFrame
        if len(Mdb) > 0:
            # These should never be NA
            int_cols = ['position', 'mm']
            for c in int_cols:
                Mdb[c] = Mdb[c].astype(int)

            # Handle bools
            bool_cols = ['consensus_SNP', 'population_SNP']
            for c in bool_cols:
                Mdb[c] = Mdb[c].astype(bool)

        self.comparison_db = pd.concat(cdbs, sort=False)
        self.mismatch_location_db = Mdb
        self.scaff2pair2mm2overlap = scaff2pair2mm2overlap

    def store_results(self):
        # Store the results in the RC
        inStrain.logUtils.log_checkpoint("Compare", "SaveResults", "start")

        s2l = {}
        for SC in self.SC_objects:
            s2l[SC.scaffold] = SC.length

        self.RCprof.store('comparisonsTable', self.comparison_db, 'pandas', 'Comparisons between the requested IS objects')
        self.RCprof.store('scaffold2length', s2l, 'dictionary', 'Scaffold to length')

        # Make the output files
        self.RCprof.generate('comparisonsTable')

        # Store scaff2pair2mm2SNPs
        if self.args.store_mismatch_locations:
            self.RCprof.store('pairwise_SNP_locations', self.mismatch_location_db, 'pandas',
                         'A dataframe of scaffold, IS pair, mm level, SNP locations')
            self.RCprof.generate('pairwise_SNP_locations')

        # Store scaff2pair2mm2cov
        if self.args.store_coverage_overlap:
            self.RCprof.store('scaff2pair2mm2cov', self.scaff2pair2mm2overlap, 'special',
                         'A dictionary of scaffold -> IS pair -> mm level -> positions with coverage overlap')

        inStrain.logUtils.log_checkpoint("Compare", "SaveResults", "end")

class ScaffoldComparison(object):
    """
    Holds the information needed to run a scaffold comparison and the results of said comparison
    """
    def __init__(self, scaffold, length):
        """
        scaffold = name of the scaffold
        length = length of scaffold
        profiles = list of ISP profiles with this scaffold (same order as names)
        names = list of names of ISP profiles with this scaffold (same order as profiles)
        """
        self.scaffold = scaffold
        self.length = int(length)

        self.profiles = []
        self.names = []

        self.SNPtables = []
        self.covTs = []

    def add_profile(self, ISP, name):
        """
        Add an ISP object with this scaffold that should be compared
        """
        self.profiles.append(ISP)
        self.names.append(name)

    def valid(self):
        """
        Return True if this scaffold has more than 2 profiles in it that can be compared
        """
        assert len(self.names) == len(set(self.names)), "You're comparing ISP profiles that have the same names. No no."
        return len(self.profiles) > 1

    def compare(self, OSC):
        """
        Return the similarity of this scaffold to the other, based on overname of names
        """
        n1 = set(self.names)
        n2 = set(OSC.names)
        return len(n1.intersection(n2)) / len(n1.union(n2))

class ScaffoldCompareGroup(object):
    """
    Holds a group of scaffolds that will be parallelized together. Ideally these scaffolds should have a high overlap
    in the profiles that they are contained in
    """
    def __init__(self, SCs):
        """
        SCs = list of scaffold comparison objects
        """
        self.ScaffoldComparisons = SCs

        # Establish the ISPs and names associated with this group of scaffolds
        ISPs = []
        names = []
        scaffolds = []
        for sc in SCs:
            scaffolds.append(sc.scaffold)
            for isp, name in zip(sc.profiles, sc.names):
                if name not in names:
                    ISPs.append(isp)
                    names.append(name)

        self.ISPs = ISPs
        self.names = names
        self.scaffolds = scaffolds

    def load_cache(self):
        """
        Load information from ISPs from disk to this object
        """
        names = self.names
        sProfiles = self.ISPs
        scaffolds_to_compare = set(self.scaffolds)

        # Make sure no duplicates in names
        assert len(names) == len(set(names)), 'Cannot have 2 of the same named IS {0}'.format(names)

        # Load data from the profiles
        name2covT = {}
        name2SNPtable = {}

        for S, name in zip(sProfiles, names):
            name2covT[name] = S.get('covT', scaffolds=scaffolds_to_compare)
            name2SNPtable[name] = S.get('cumulative_snv_table').rename(
                columns={'conBase': 'con_base', 'refBase': 'ref_base', 'varBase': 'var_base',
                         'baseCoverage': 'position_coverage'})

        # Attach this information to ScaffoldComparison objects
        for SC in self.ScaffoldComparisons:
            for name in SC.names:
                SC.SNPtables.append(inStrain.compare_utils.subset_SNP_table(name2SNPtable[name], SC.scaffold))
                SC.covTs.append(name2covT[name][SC.scaffold])

        # Set up a command and result queue
        ctx = multiprocessing.get_context('spawn')
        self.cmd_queue = ctx.Queue()
        self.result_queue = ctx.Queue()
        for SC in self.ScaffoldComparisons:
            self.cmd_queue.put(SC)

    def purge_cache(self):
        """
        Purge information from ISPs from RAM
        """
        for SC in self.ScaffoldComparisons:
            del SC.SNPtables
            del SC.covTs

    def __str__(self):
        """
        String representation
        """
        string = f"This group has {len(self.ScaffoldComparisons)} scaffolds accross {len(self.names)} ISP objects"
        return string

def group_Scaffold_objects(valid_SCs, group_length):
    SC_groups = []
    current_group = []
    current_size = 0
    for s in valid_SCs:
        current_group.append(s)
        current_size += s.length
        if current_size >= group_length:
            SC_groups.append(ScaffoldCompareGroup(current_group))
            current_group = []
            current_size = 0

    if current_size > 0:
        SC_groups.append(ScaffoldCompareGroup(current_group))

    return SC_groups

def make_scaffoldcomparison_objects(inputs, scaffolds_to_compare):
    inStrain.logUtils.log_checkpoint("Compare", "CreateScaffoldComparisonObjects", "start")

    # Get the stuff to return
    scaffold2SC = {}

    # Keep this cached; necessary for genome-level operations
    scaffold2length = {}

    # Go through the profiles
    for profile_loc in inputs:
        if not os.path.exists(profile_loc):
            logging.error("IS {0} does not exist! Skipping".format(profile_loc))
            continue

        logging.info("Loading {0}".format(profile_loc))
        ISP = inStrain.SNVprofile.SNVprofile(profile_loc)
        scaffolds = list(ISP._get_covt_keys())
        name = os.path.basename(ISP.get('bam_loc'))

        if scaffolds_to_compare is not None:
            scaffolds = list(set(scaffolds).intersection(scaffolds_to_compare))

        # Update scaffold to length
        s2l = ISP.get('scaffold2length')
        for scaff, l in s2l.items():
            if scaff in scaffold2length:
                assert l == scaffold2length[scaff]
            else:
                scaffold2length[scaff] = l

        # Create ScaffoldComparison objects
        for scaff in scaffolds:
            if scaff not in scaffold2SC:
                SC = ScaffoldComparison(scaff, scaffold2length[scaff])
                scaffold2SC[scaff] = SC
            scaffold2SC[scaff].add_profile(ISP, name)

    # Figure out which scaffolds to compare
    valid_SCs = [SC for scaff, SC in scaffold2SC.items() if SC.valid()]
    logging.info(f"{len(valid_SCs)} of {len(scaffold2SC.keys())} scaffolds are in at least 2 samples")
    assert len(valid_SCs) > 0, "No scaffolds are shared among the IS objects"

    inStrain.logUtils.log_checkpoint("Compare", "CreateScaffoldComparisonObjects", "end")

    return valid_SCs, scaffold2length

def calc_scaff_sim_matrix(valid_SCs):
    """
    Perform pairwise comparison of scaffolds based on the overlap of profiles
    """
    table = defaultdict(list)
    for i, sc1 in enumerate(valid_SCs):
        for j, sc2 in enumerate(valid_SCs):
            if i <= j:
                continue
            table['sc1'].append(sc1)
            table['sc2'].append(sc2)
            table['similarity'].append(sc1.compare(sc2))
    return pd.DataFrame(table)

def establish_SC_groups(valid_SCs, Cdb, group_length):
    """
    Based on the distance matrix, put valid_SCs into a number of groups with a max of group_length
    """
    # Just do it simple; would be better to use the distance matrix, but oh well
    SC_groups = simple_grouping(valid_SCs, group_length)

    return SC_groups