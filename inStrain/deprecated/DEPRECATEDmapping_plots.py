"""
I tried to go down a seaborn-inspired structure here, but gave up on it. Leaving here in case want to go down
this road again in the future.

- MO 2/3/21
"""


import logging
import traceback

import matplotlib
matplotlib.use('Agg')
import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

import inStrain.plotting
import inStrain.plotting.utilities

def BreadthCurvePlot(IS, **kwargs):
    """
    Create a breadth curve plot (plot number 1)
    """
    plotter = _BreadthCurvePlotter(IS, **kwargs)

    logging.info("Plotting plot 1")
    plotter.plot()

class _BreadthCurvePlotter:
    """
    Draw the breadth curve plot
    """
    def plot(self):
        logging.info("Plotting plot 1")

        # Establish multipage .pdf
        pp = PdfPages(self.plot_dir + self.name)

        # Make the plots
        for genome, mdb in self.GWdb.groupby('genome'):

            # See if you should plot this genome
            if not inStrain.plotting.utilities.plot_genome(genome, self.IS, **self.kwargs):
                continue

    def __init__(self, IS, **kwargs):
        """
        Initialize the plotting object
        """
        # Load simple kwargs
        self.debug = kwargs.get('debug', False)
        self.plot_dir = kwargs.get('plot_dir', None)
        self.name = kwargs.get('name', 'CoverageAndBreadth_vs_readMismatch.pdf')
        self.kwrags = kwargs
        self.IS = IS

        # Load data
        try:
            # Load GWdb cache
            if kwargs.get('GWdb', None) is not None:
                GWdb = kwargs.get('GWdb')
            else:
                GWdb = None
                print('Make GWdb you fool')
            assert len(GWdb) > 0

            # Parse GWdb a bit
            readLen = int(IS.get_read_length())
            GWdb['read_length'] = readLen
            GWdb['mm'] = GWdb['mm'].astype(int)
            GWdb['ANI_level'] = [(readLen - mm) / readLen for mm in GWdb['mm']]

        except:
            logging.error(
                "Skipping plot 1 - you don't have all required information.")
            traceback.print_exc()
            return

