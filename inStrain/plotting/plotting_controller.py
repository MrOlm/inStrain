import os
import sys
import logging
import traceback

import inStrain
import inStrain.SNVprofile
import inStrain.plotting
import inStrain.plotting.utilities

from inStrain.plotting.mapping_plots import mm_plot_from_IS, ANI_dist_plot_from_IS, read_filtering_from_IS
from inStrain.plotting.positional_plots import genome_plot_from_IS, scaffold_inspection_from_IS
from inStrain.plotting.SNV_plots import allele_freq_plot_from_IS
from inStrain.plotting.gene_plots import gene_histogram_from_IS
from inStrain.plotting.linkage_plots import linkage_decay_from_IS, linkage_decay_type_from_IS
from inStrain.plotting.compare_plots import dendrograms_from_RC

class PlottingController(object):
    """
    Handle the logic of profiling
    """
    PLOTTING_FUNCTIONS = [
        [1, mm_plot_from_IS],
        [2, genome_plot_from_IS],
        [3, ANI_dist_plot_from_IS],
        [4, allele_freq_plot_from_IS],
        [5, linkage_decay_from_IS],
        [6, read_filtering_from_IS],
        [7, scaffold_inspection_from_IS],
        [8, linkage_decay_type_from_IS],
        [9, gene_histogram_from_IS],
        [10, dendrograms_from_RC]
    ]

    def __init__(self, args):
        """
        Accepts args directly from argparse
        """
        self.ori_args = args
        self.args = args

        self.validate_input()

    def main(self):
        """
        Profile the bam with inStrain using the options in kwargs
        """
        if not self.validated:
            return

        # Establish kwargs, which is kind of like the cacheing of this module
        self.establish_kwargs()

        # Run the plots
        self.run_plots()

    def run_plots(self):
        """
        Call to create all the plots
        """
        for plf in PlottingController.PLOTTING_FUNCTIONS:
            number, func = plf
            if str(number) in self.to_plot:
                try:
                    func(self.IS, plot_dir=self.plot_dir, **self.kwargs)
                except BaseException as e:
                    logging.error(f'Failed to make plot #{number}: {str(e)}')
                    if self.debug:
                        traceback.print_exc()
                        logging.debug(traceback.format_exc())


    def establish_kwargs(self):
        """
        Load data from the IS file that will be based to plots through the kwargs
        """
        # Remove IS from the kwargs
        kwargs = vars(self.args)
        kwargs.pop('IS')

        # Cache needed data
        if self.IS_TYPE == 'IS':
            try:
                kwargs['GWdb'] = self.IS.get('genome_level_info')
            except:
                logging.error(
                    "Cannot cache scaffold info - you don't have all required information. You need to run inStrain genome_wide first")
                if self.debug:
                    traceback.print_exc()

        if (('2' in self.to_plot) | ('7' in self.to_plot)):
            kwargs['covT'] = self.IS.get('covT')
            kwargs['clonT'] = self.IS.get('clonT')
            kwargs['raw_linkage_table'] = self.IS.get('raw_linkage_table')
            kwargs['cumulative_snv_table'] = self.IS.get('cumulative_snv_table')

        self.kwargs = kwargs

    def validate_input(self):
        """
        Do a bunch of parsing of the input arguments
        """
        args = self.args
        kwargs = vars(args)

        # Get the IS object
        assert os.path.exists(args.IS)
        self.IS = inStrain.SNVprofile.SNVprofile(args.IS)

        # Set up the logger
        log_loc = self.IS.get_location('log') + 'log.log'
        inStrain.controller.setup_logger(log_loc)
        logging.getLogger('matplotlib.font_manager').disabled = True

        # Set up the .stb
        stb = self.IS.get('scaffold2bin')
        if stb is None:
            logging.error("This IS object does not have an .stb file; cant use it to make plots")
            self.validated = False

        # Set debug
        self.debug = args.debug

        # Figure out if this is an RC or a IS
        bl = self.IS.get('bam_loc')
        if bl is None:
            self.IS_TYPE = 'RC'
            self.options = ['10']
        else:
            self.IS_TYPE = 'IS'
            self.options = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

        # Get the plot directory and basename
        self.plot_dir = self.IS.get_location('figures') + os.path.basename(self.IS.get('location')) + '_'

        # Figure out which plots to make
        self.to_plot = _parse_plot_options(self.options, kwargs.get('plots', None))
        logging.info("making plots {0}".format(', '.join(self.to_plot)))
        self.validated = True

def _parse_plot_options(options, args):
    '''
    Read user input and figure out a list of plots to make

    Args:
        options: list of possible plots to make (default [1-6])
        args: the command line passed in

    Returns:
        list: list of ints in the args
    '''
    to_plot = []

    if args[0] in ['all', 'a']:
        to_plot += options

    elif args == None:
        logging.error("No plots given!")
        return []

    else:
        for arg in args:
            if arg in options:
                to_plot.append(arg)

    return to_plot

