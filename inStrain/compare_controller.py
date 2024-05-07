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
import traceback
from collections import defaultdict
from tqdm import tqdm
import time

import inStrain.readComparer
import inStrain.compare_utils
import inStrain.genomeUtilities
import inStrain.controller
import inStrain.profile
import inStrain.profile.samtools_ops

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
        message = """\
***************************************************
    ..:: inStrain compare Step 1. Load data ::..
***************************************************
                """
        logging.info(message)
        self.create_scaffoldcomparison_objects()

        # Actually do the processing
        message = """\
***************************************************
..:: inStrain compare Step 2. Run comparisons ::..
***************************************************
                """
        logging.info(message)
        self.run_comparisons()

        # Do auxillary processing if needed
        message = """\
***************************************************
..:: inStrain compare Step 3. Auxiliary processing ::..
***************************************************
                """
        logging.info(message)
        self.run_auxillary_processing()

        # Store the results
        message = """\
***************************************************
..:: inStrain compare Step 4. Store results ::..
***************************************************
                """
        logging.info(message)
        self.store_results()
        self.write_final_message()

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
                # Load it
                stb = inStrain.genomeUtilities.load_scaff2bin(args.stb)
                self.stb = stb

                # Only compare scaffolds in the .stb
                if self.scaffolds_to_compare is not None:
                    self.scaffolds_to_compare = self.scaffolds_to_compare.union(set(stb.keys()))
                else:
                    self.scaffolds_to_compare = set(stb.keys())

        # Load database mode results
        if self.args.database_mode is True:
            inStrain.logUtils.log_checkpoint("Compare", "LoadDatabaseMode", "start")

            bin2scaffolds = defaultdict(list)
            for s, b in stb.items():
                bin2scaffolds[b].append(s)

            input2scaffolds = {}
            for input in self.inputs:
                scaffolds = inStrain.compare_utils.find_relevant_scaffolds(input, bin2scaffolds, self.kwargs)
                input2scaffolds[input] = scaffolds

            self.input2scaffolds = input2scaffolds
            inStrain.logUtils.log_checkpoint("Compare", "LoadDatabaseMode", "end")

        # Handle single genome mode
        if self.args.genome is not None:
            bin2scaffolds = defaultdict(list)
            for s, b in stb.items():
                bin2scaffolds[b].append(s)
            if self.args.genome in bin2scaffolds:
                self.scaffolds_to_compare = self.scaffolds_to_compare.intersection(set(bin2scaffolds[self.args.genome]))
            else:
                logging.error(f'genome {self.args.genome} is not in the provided .stb file')
                raise Exception

        # Load the .bam files
        self.name2bam = gen_name_2_bam(self.inputs, args.bams, self.kwargs.get('processes', 1))

        i = 0
        if self.name2bam != {}:
            self.run_pooling = True
            for name, bam in self.name2bam.items():
                input = self.inputs[i]
                i += 1

                if name != os.path.basename(bam):
                    logging.error(f"WARNING - .bam file {os.path.basename(bam)} is being used for input {input}"\
                                  f" even though the .bam used in creation of that input is {name}. Change order of "\
                                  f"--bams input to fix this if it's not correct")
        else:
            self.run_pooling = False

    def create_scaffoldcomparison_objects(self):
        """
        Figure out which scaffolds need to be compared
        """
        # Load scaffolds to compare
        if self.scaffolds_to_compare is not None:
            scaffolds_to_compare = set(self.scaffolds_to_compare)
        else:
            scaffolds_to_compare = None

        # Load input2scaffolds
        if self.args.database_mode is True:
            i2s = self.input2scaffolds
        else:
            i2s = None

        # Create comparison objects
        valid_SCs, scaffold2length = make_scaffoldcomparison_objects(self.inputs, scaffolds_to_compare,
                                                                     input2scaffolds=i2s)

        self.SC_objects = valid_SCs
        self.scaffold2length = scaffold2length

        # Establish ScaffoldComparison groups
        group_length = self.kwargs.get('group_length', 10000000)

        #Cdb = calc_scaff_sim_matrix(valid_SCs)
        #SC_groups = establish_SC_groups(valid_SCs, Cdb, group_length)
        SC_groups = group_Scaffold_objects(valid_SCs, group_length)

        self.scaffold_comparison_groups = SC_groups

    def run_comparisons(self):
        """
        Run each group one at a time, storing the results in the meantime
        """
        if self.kwargs.get('skip_popANI', False):
            return

        cdbs = [] # Store scaffold comparison information
        mdbs = [] # Store mismatch locations
        pair2mm2covOverlaps = [] # Store coverage overlap locations
        order = [] # Store order of scaffolds

        groups = len(self.scaffold_comparison_groups)

        inStrain.logUtils.log_checkpoint("Compare", "multiprocessing", "start")
        for i, SCgroup in enumerate(self.scaffold_comparison_groups):
            logging.info(f'Running group {i+1} of {groups}')

            SCgroup.load_cache()
            results = inStrain.compare_utils.run_compare_multiprocessing(SCgroup.cmd_queue, SCgroup.result_queue,
                                                                         self.null_model, num_to_run=len(SCgroup.scaffolds),
                                                                         **self.kwargs)
            for result in results:
                if result is not None:
                    # This means you're not doing any pooling
                    assert len(result) == 4
                    Cdb, Mdb, pair2mm2covOverlap, scaffold = result
                    for item, lis in zip([Cdb, Mdb, pair2mm2covOverlap, scaffold], [cdbs, mdbs, pair2mm2covOverlaps, order]):
                        lis.append(item)

            SCgroup.purge_cache()

        # Process main results
        inStrain.logUtils.log_checkpoint("Compare", "multiprocessing", "end")
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

    # def process_pooling_results(self, DSTdbs, PMdbs, order):
    #     """
    #     Merge and store pooling results
    #     """
    #     # Create a "scaffold" categorical
    #     # https://pandas.pydata.org/pandas-docs/stable/user_guide/categorical.html
    #     from pandas.api.types import CategoricalDtype
    #     scaff_cat = CategoricalDtype(order)
    #
    #     # Add this category to all the dataframes
    #     for d, s in zip(DSTdbs, order):
    #         d['scaffold'] = s
    #         d['scaffold'] = d['scaffold'].astype(scaff_cat)
    #
    #     # Add this category to all the dataframes
    #     for d, s in zip(PMdbs, order):
    #         d['scaffold'] = s
    #         d['scaffold'] = d['scaffold'].astype(scaff_cat)
    #
    #     self.DSTdb = pd.concat(DSTdbs)
    #     self.PMdb = pd.concat(PMdbs).astype(
    #                 {'A':int, 'C':int, 'T':int, 'G':int, 'depth':int,
    #                  'sample_detections':int, 'DivergentSite_count':int, 'SNS_count':int, 'SNV_count':int,
    #                  'con_SNV_count':int, 'pop_SNV_count':int, 'sample_5x_detections':int})

    def run_auxillary_processing(self):
        """
        Handle .stb files, database mode, and more
        """
        # Run pooling
        if self.run_pooling:
            PM = inStrain.polymorpher.PoolController(self.SC_objects, self.name2bam)
            PM.main()
            self.DSTdb = PM.DSTdb
            self.PMdb = PM.PMdb

        # Calculate s2l (scaffold 2 length)
        s2l = {}
        for SC in self.SC_objects:
            s2l[SC.scaffold] = SC.length
        self.s2l = s2l

        # Calculate genome-level results
        if hasattr(self, 'stb') and not self.kwargs.get('skip_popANI', False):
            # Calculate bin 2 length
            # Make bin to length
            b2l = {}
            for scaffold, bin in self.stb.items():
                if bin not in b2l:
                    b2l[bin] = 0

                if scaffold in s2l:
                    b2l[bin] += s2l[scaffold]
                else:
                    logging.debug(
                        "FAILURE StbError {0} {1} no_length will not be considered as part of the genome".format(
                            scaffold, bin))
            self.bin2length = b2l

            gdb = inStrain.genomeUtilities._add_stb(self.comparison_db, self.stb)
            Gdb = inStrain.genomeUtilities._genome_wide_readComparer(gdb, self.stb, b2l, **self.kwargs)
            self.genomelevel_compare = Gdb

        # Cluster genomes
        if hasattr(self, 'genomelevel_compare') and not self.kwargs.get('skip_popANI', False):
            self.run_genome_clustering()

    def store_results(self):
        # Store the results in the RC
        inStrain.logUtils.log_checkpoint("Compare", "SaveResults", "start")

        if hasattr(self, 'comparison_db'):
            self.RCprof.store('comparisonsTable', self.comparison_db, 'pandas', 'Comparisons between the requested IS objects')

        if hasattr(self, 's2l'):
            self.RCprof.store('scaffold2length', self.s2l, 'dictionary', 'Scaffold to length')

        force_compress = vars(self.args).get('force_compress', False)

        # Store auxillary things
        if hasattr(self, 'bin2length'):
            self.RCprof.store('bin2length', self.bin2length, 'dictionary', 'Dictionary of bin 2 total length')
        if hasattr(self, 'genomelevel_compare'):
            # ... this is jankey, but I guess it's what we do
            out_base = self.RCprof.get_location('output') + os.path.basename(self.RCprof.get('location')) + '_'
            name = out_base + 'genomeWide_compare.tsv'
            if force_compress:
                name += '.gz'
            self.genomelevel_compare.to_csv(name, index=False, sep='\t')
        if hasattr(self, 'Cdb'):
            # ... this is jankey, but I guess it's what we do
            out_base = self.RCprof.get_location('output') + os.path.basename(self.RCprof.get('location')) + '_'
            name = out_base + 'strain_clusters.tsv'
            if force_compress:
                name += '.gz'
            self.Cdb.to_csv(name, index=False, sep='\t')
        if hasattr(self, 'stb'):
            self.RCprof.store('scaffold2bin', self.stb, 'dictionary', 'Dictionary of scaffold 2 bin')

        # Make the output files
        if not self.kwargs.get('skip_popANI', False):
            self.RCprof.generate('comparisonsTable', **vars(self.args))

        # Store scaff2pair2mm2SNPs
        if self.args.store_mismatch_locations:
            self.RCprof.store('pairwise_SNP_locations', self.mismatch_location_db, 'pandas',
                         'A dataframe of scaffold, IS pair, mm level, SNP locations')
            self.RCprof.generate('pairwise_SNP_locations', **vars(self.args))

        # Store scaff2pair2mm2cov
        if self.args.store_coverage_overlap:
            self.RCprof.store('scaff2pair2mm2cov', self.scaff2pair2mm2overlap, 'special',
                         'A dictionary of scaffold -> IS pair -> mm level -> positions with coverage overlap')

        # Store pooling
        if hasattr(self, 'DSTdb'):
            self.RCprof.store('DSTdb', self.DSTdb, 'pickle',
                             'Dense SNV table (from pooling)')

        if hasattr(self, 'PMdb'):
            self.RCprof.store('PMdb', self.PMdb, 'pickle',
                             'Informative SNV table (from pooling)')

        if self.run_pooling and hasattr(self, 'PMdb'):
            self.RCprof.generate('pooled_SNV_info', **vars(self.args))
            self.RCprof.generate('pooled_SNV_data', **vars(self.args))

        inStrain.logUtils.log_checkpoint("Compare", "SaveResults", "end")

        # Make plots
        if not self.kwargs.get('skip_popANI', False):
            if hasattr(self, 'genomelevel_compare'):
                self.make_plots()

    def make_plots(self):
        '''
        Call plotting function from "profile" module
        '''
        args = self.args
        args.IS = self.RCprof.location

        if not args.skip_plot_generation:
            inStrain.logUtils.log_checkpoint("Compare", "making_plots", "start")
            args.plots = ['10']
            inStrain.controller.Controller().plot_operation(args)
            inStrain.logUtils.log_checkpoint("Compare", "making_plots", "end")
        else:
            logging.info('Nevermind! You chose to skip making plots')

    def write_final_message(self):
        Sprofile = self.RCprof
        message = """\
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

..:: inStrain compare finished ::..

Output tables........ {0}
Figures.............. {1}
Logging.............. {2}

See documentation for output descriptions - https://instrain.readthedocs.io/en/latest/

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        """.format(Sprofile.get_location('output'), \
                   Sprofile.get_location('figures'),
                   Sprofile.get_location('log'))
        logging.info(message)

    def run_genome_clustering(self):
        """
        Cluster genomes dRep style
        """
        kwargs = self.kwargs
        Mdb = self.genomelevel_compare
        Mdb = Mdb.sort_values(['genome', 'name1', 'name2'])
        if len(Mdb) == 0:
            return

        try:
            Cdb = inStrain.compare_utils.cluster_genome_strains(Mdb, kwargs)
            self.Cdb = Cdb
        except:
            logging.error('Could not cluster genomes; heres a traceback:')
            traceback.print_exc()

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
        assert len(self.names) == len(set(self.names)), "InStrain profiles that have the same name cannot be compared. See https://github.com/MrOlm/inStrain/issues/79 and https://github.com/MrOlm/inStrain/issues/97"
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

    def load_cache(self, for_pooling=False):
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
        name2SNPtable_hash = {}
        name2Rdic = {}

        for S, name in zip(sProfiles, names):
            logging.debug('Loading %s', name)
            start = time.time()
            if not for_pooling:
                name2covT[name] = S.get('covT', scaffolds=scaffolds_to_compare)
            logging.debug('covT %f', time.time() - start)
            start = time.time()
            name2SNPtable[name] = S.get('cumulative_snv_table').rename(
                columns={'conBase': 'con_base', 'refBase': 'ref_base', 'varBase': 'var_base',
                         'baseCoverage': 'position_coverage'})
            print(name)
            name2SNPtable_hash[name] = inStrain.compare_utils.hash_SNP_table(name2SNPtable[name])
            logging.debug('cumulative_snv_table %f', time.time() - start)
            if for_pooling:
                start = time.time()
                name2Rdic[name] = S.get("Rdic")
                logging.debug('Rdic %f', time.time() - start)

        if for_pooling:
            self.name2SNPtable = name2SNPtable
            self.name2Rdic = name2Rdic

        if not for_pooling:
            # Attach this information to ScaffoldComparison objects
            for SC in self.ScaffoldComparisons:
                for name in SC.names:
                    # SC.SNPtables.append(inStrain.compare_utils.subset_SNP_table(name2SNPtable[name], SC.scaffold))
                    SC.SNPtables.append(name2SNPtable_hash[name].get(SC.scaffold, pd.DataFrame()))
                    SC.covTs.append(name2covT[name][SC.scaffold])

                    if for_pooling:
                        SC.name2Rdic = {n:r[SC.scaffold] for n, r in name2Rdic.items()}

            # Set up a command and result queue
            ctx = multiprocessing.get_context('spawn')
            self.cmd_queue = ctx.Queue()
            self.result_queue = ctx.Queue()
            for SC in self.ScaffoldComparisons:
                self.cmd_queue.put(SC)

        logging.debug('Load cache done.')

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

def make_scaffoldcomparison_objects(inputs, scaffolds_to_compare, input2scaffolds=None):
    inStrain.logUtils.log_checkpoint("Compare", "CreateScaffoldComparisonObjects", "start")

    # Get the stuff to return
    scaffold2SC = {}

    # Keep this cached; necessary for genome-level operations
    scaffold2length = {}

    # Go through the profiles
    for profile_loc in tqdm(inputs, desc='Loading Profiles into RAM'):
        if not os.path.exists(profile_loc):
            logging.error("IS {0} does not exist! Skipping".format(profile_loc))
            continue

        logging.debug("Loading {0}".format(profile_loc))
        ISP = inStrain.SNVprofile.SNVprofile(profile_loc)
        scaffolds = list(ISP._get_covt_keys())
        name = os.path.basename(ISP.get('bam_loc'))

        if input2scaffolds is not None:
            scaffolds = list(set(scaffolds).intersection(input2scaffolds[profile_loc]))
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

def gen_name_2_bam(inputs, bams, processes=1):
    name2bam = {}

    # See if you have a list of bam objects
    if type(bams) == type([]):
        if len(bams) == len(inputs):
            for profile_loc, bam in zip(inputs, bams):
                ISP = inStrain.SNVprofile.SNVprofile(profile_loc)
                name = os.path.basename(ISP.get('bam_loc'))

                name2bam[name] = inStrain.profile.samtools_ops.prepare_bam_fie(bam, processes)

    return name2bam
