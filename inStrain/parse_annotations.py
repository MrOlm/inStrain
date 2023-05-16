import os
import copy
import logging
import pandas as pd

from collections import defaultdict
from tqdm import tqdm

import inStrain
import inStrain.SNVprofile

class PAController(object):
    '''
    Main controller of the parse_annotations command
    '''
    def __init__(self, args):
        '''
        Set all of the command line arguments in the "args" attribute
        '''
        self.ori_args = copy.deepcopy(args)
        self.kwargs = vars(args)

    def main(self):
        '''
        The main method when run on the command line
        '''
        # Validate arguments and such
        self.parse_arguments()

        # Load the annotation table
        message = """\
***************************************************
    ..:: Step 1. Load tables ::..
***************************************************
                """
        logging.info(message)
        self.load_annotation_table_wrapper()
        self.load_gene_tables_wrapper()

        # Run the calculations
        message = """\
***************************************************
..:: Step 2. Run calculations ::..
***************************************************
                """
        logging.info(message)
        self.gene_calculations_wrapper()

        # Store the results
        message = """\
***************************************************
..:: Step 3. Store results ::..
***************************************************
                """
        logging.info(message)
        self.store_results()
        self.write_final_message()
#
    def parse_arguments(self):
        """
        Parse the arguments and add them this object as attributes
        """
        args = self.kwargs

        # Set up output object and log
        outbase = args.get('output')
        RCprof = inStrain.SNVprofile.SNVprofile(outbase)
        log_loc = RCprof.get_location('log') + 'log.log'
        inStrain.controller.setup_logger(log_loc)
        self.OD = RCprof

        # Set up list of input IS profiles
        inputs = args.get('input')
        assert len(inputs) >= 1, "You need to have at least one input .IS file"
        self.input_ISs = [inStrain.SNVprofile.SNVprofile(i) for i in inputs]

    def load_annotation_table_wrapper(self):
        """
        Wrapper for loading the annotation table
        """
        args = self.kwargs
        gene2annos = load_annotation_table2(args.get('annotations'))
        self.gene2anno = gene2annos

    def load_gene_tables_wrapper(self):
        """
        Wrapper for loading the gene tables
        """
        args = self.kwargs

        gdbs = []
        names = []
        for IS in self.input_ISs:
            gdb = IS.load_output('gene_info')
            gdb['annos'] = gdb['gene'].map(self.gene2anno)
            name = os.path.basename(IS.get('bam_loc'))

            if args.get('min_genome_breadth') > 0:
                genome_db = IS.load_output('genome_info')
                stb = IS.get('scaffold2bin')
                if len(genome_db) > 0:
                    genomes = set(genome_db[genome_db['breadth'] >= args.get('min_genome_breadth')]['genome'])
                else:
                    logging.error("No genomes detected in this sample!")
                    genomes = []
            else:
                genomes = None
                stb = None

            gdb = filter_gene_table(gdb, genomes, stb, **self.kwargs)
            gdbs.append(gdb)
            names.append(name)

        self.gene_tables = gdbs
        self.names = names

    def gene_calculations_wrapper(self):
        # Calculate summary stats
        sdb = calculate_gene_sum_stats(self.gene_tables, self.names, **self.kwargs)
        self.totals_db = sdb

        # Calculate annotation counts
        s2a2g2vals = calculate_annotation_counts2(self.gene_tables, self.names, **self.kwargs)

        # Create tables
        metric2table = create_annotation_tables2(sdb, s2a2g2vals, **self.kwargs)
        self.metric2tables = metric2table

    def store_results(self):
        # Store raw data
        if self.kwargs.get('store_rawdata', False):
            self.OD.store('gene2anno', self.gene2anno, 'dictionary', 'Dictionary of genes 2 annotations')

        # Store output tables
        outloc = self.OD.get_location('output')

        self.totals_db.to_csv(outloc + 'SampleAnnotationTotals.csv', index=False)

        m2n = {'long_data':'LongFormData.csv'}
        for metric, table in self.metric2tables.items():
            if metric in m2n:
                name = m2n[metric]
            else:
                name = 'ParsedGeneAnno_' + metric + '.csv'
            table.to_csv(outloc + name, index=False)

    def write_final_message(self):
        Sprofile = self.OD
        message = """\
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

..:: inStrain parse_annotations finished ::..

Output tables........ {0}
Logging.............. {1}

See documentation for output descriptions - https://instrain.readthedocs.io/en/latest/

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                """.format(Sprofile.get_location('output'), \
                           Sprofile.get_location('log'))
        logging.info(message)

def load_annotation_table2(locs):
    """
    Load the annotation table

    Input can either be a list of file locations or just 1 file location
    """
    gene2annos = {}

    if type(locs) != type([]):
        locs = [locs]

    for loc in tqdm(locs, total=len(locs), desc='Loading sample gene data'):
        first = True
        with open(loc) as o:
            for line in o:
                lw = [l.strip() for l in line.split(',')]
                if first:
                    first = False
                    for c in ['gene', 'anno']:
                        if c not in list(lw):
                            m = f"The required column {c} is not in your annotation table {loc}!"
                            print(m)
                            logging.error(m)
                            raise Exception(m)
                    continue

                if len(lw) != 2:
                    m = f"The line {line} cannot be parsed!"
                    print(m)
                    logging.error(m)
                    continue

                gene = lw[0]
                anno = lw[1]
                if gene not in gene2annos:
                    gene2annos[gene] = set()
                gene2annos[gene].add(anno)

    return gene2annos

def filter_gene_table(gdb, genomes=None, stb=None, min_gene_breadth=0.5, **kwargs):
    """
    Filter the gene table according to the kwargs
    """
    db = gdb[gdb['coverage'] > 0]
    db = db[db['breadth'] >= min_gene_breadth]

    if genomes is not None:
        db['genome'] = db['scaffold'].map(stb)
        db = db[db['genome'].isin(genomes)]

    return db

def calculate_gene_sum_stats(gdbs, names, **kwargs):
    """
    Calculate summary stats
    """
    table = defaultdict(list)
    for db, name in zip(gdbs, names):
        db['gene_length'] = abs(db['end'] - db['start']) + 1
        db['mb'] = db['gene_length'] * db['coverage']
        db['mb'] = db['mb'].astype(int)

        table['sample'].append(name)
        table['detected_genes'].append(len(db))

        if 'genome' in set(db.columns):
            table['detected_genomes'].append(len(db['genome'].unique()))

        table['bases_mapped_to_genes'].append(int(db['mb'].sum()))
        table['detected_annotations'].append(sum([len(a) for a in db['annos'] if a == a]))
        table['detected_genes_with_anno'].append(len(db[~db['annos'].isna()]))

    return pd.DataFrame(table)

# def calculate_annotation_counts(gdbs, names, **kwargs):
#     """
#     Create sample2accession2values
#     """
#     s2a2vals = {}
#     for db, name in tqdm(zip(gdbs, names), total=len(names), desc='Calculating anno metrics'):
#         # Create "accession to values"
#         # values = [genomes, # genes, # bases]
#         a2vals = {}
#         for i, row in db[~db['annos'].isna()].iterrows():
#             kos = row['annos']
#             if 'genome' in row:
#                 genome = row['genome']
#             else:
#                 genome = None
#             for k in kos:
#                 if k in a2vals:
#                     a2vals[k][0] = a2vals[k][0].union(set([genome]))
#                     a2vals[k][1] += 1
#                     a2vals[k][2] += row['mb']
#                 else:
#                     a2vals[k] = [set([genome]), 1, row['mb']]
#         s2a2vals[name] = a2vals
#     return s2a2vals

def calculate_annotation_counts2(gdbs, names, **kwargs):
    """
    Create sample2accession2genome2values
    """
    s2a2g2vals = {}
    for db, name in tqdm(zip(gdbs, names), total=len(names), desc='Calculating anno metrics'):
        # Create "accession to genome to values"
        # values = [genomes, # genes, # bases]
        a2g2vals = {}
        for i, row in db[~db['annos'].isna()].iterrows():
            kos = row['annos']
            if 'genome' in row:
                g = row['genome']
            else:
                g = None
            for k in kos:

                if k not in a2g2vals:
                    a2g2vals[k] = {}

                if g in a2g2vals[k]:
                    a2g2vals[k][g][0] = a2g2vals[k][g][0].union(set([g]))
                    a2g2vals[k][g][1] += 1
                    a2g2vals[k][g][2] += row['mb']
                else:
                    a2g2vals[k][g] = [set([g]), 1, row['mb']]
        s2a2g2vals[name] = a2g2vals
    return s2a2g2vals

# def create_annotation_tables(sdb, s2a2vals, **kwargs):
#     """
#     Create the actual tables that can be used for stats and whatnot
#     """
#     if 'detected_genomes' in set(sdb.columns):
#         METRICS = ['genes', 'bases', 'genomes']
#     else:
#         METRICS = ['genes', 'bases']
#
#     metric2table = {}
#     for metric in METRICS:
#         metric2table[metric] = defaultdict(list)
#
#     TOTAL_KOS = set()
#     for s, a2vals in s2a2vals.items():
#         TOTAL_KOS = TOTAL_KOS.union(set(a2vals.keys()))
#
#     for sample, a2vals in s2a2vals.items():
#         for metric in METRICS:
#             metric2table[metric]['sample'].append(sample)
#
#         for KO in sorted(list(TOTAL_KOS)):
#             if KO in a2vals:
#                 genomes = len(a2vals[KO][0])
#                 genes = a2vals[KO][1]
#                 bases = a2vals[KO][2]
#             else:
#                 genomes = 0
#                 genes = 0
#                 bases = 0
#
#             for metric in METRICS:
#                 if metric == 'genes':
#                     val = genes
#                 elif metric == 'bases':
#                     val = bases
#                 elif metric == 'genomes':
#                     val = genomes
#                 else:
#                     assert False
#                 metric2table[metric][KO].append(val)
#
#     for metric in METRICS:
#         metric2table[metric] = pd.DataFrame(metric2table[metric])
#
#     # Calculate the percentage metrics
#     for metric in METRICS:
#         if metric == 'genes':
#             s2norm = sdb.set_index('sample')['detected_genes'].to_dict()
#         elif metric == 'bases':
#             s2norm = sdb.set_index('sample')['bases_mapped_to_genes'].to_dict()
#         elif metric == 'genomes':
#             s2norm = sdb.set_index('sample')['detected_genomes'].to_dict()
#         else:
#             assert False
#
#         db = metric2table[metric].copy()
#         for d in TOTAL_KOS:
#             db[d] = [x / s2norm[s] for x, s in zip(db[d], db['sample'])]
#
#         metric2table[metric + '_fraction'] = db
#
#     # Calculate the long-form table
#     table = defaultdict(list)
#     for sample, a2vals in s2a2vals.items():
#         for a, vals in a2vals.items():
#             table['sample'].append(sample)
#             table['anno'].append(a)
#             for name, thing in zip(['genomes', 'genes', 'bases'], vals):
#                 table[name].append(thing)
#     metric2table["long_data"] = pd.DataFrame(table)
#
#     return metric2table

def create_annotation_tables2(sdb, s2a2g2vals, **kwargs):
    """
    Create the actual tables that can be used for stats and whatnot
    """
    if 'detected_genomes' in set(sdb.columns):
        METRICS = ['genes', 'bases', 'genomes']
    else:
        METRICS = ['genes', 'bases']

    metric2table = {}
    for metric in METRICS:
        metric2table[metric] = defaultdict(list)

    TOTAL_KOS = set()
    for s, a2g2vals in s2a2g2vals.items():
        TOTAL_KOS = TOTAL_KOS.union(set(a2g2vals.keys()))

    for sample, a2g2vals in s2a2g2vals.items():
        for metric in METRICS:
            metric2table[metric]['sample'].append(sample)

        for KO in sorted(list(TOTAL_KOS)):
            if KO in a2g2vals:
                genomes = len(a2g2vals[KO].keys())
                genes = sum([a2vals[1] for g, a2vals in a2g2vals[KO].items()])
                bases = sum([a2vals[2] for g, a2vals in a2g2vals[KO].items()])
            else:
                genomes = 0
                genes = 0
                bases = 0

            for metric in METRICS:
                if metric == 'genes':
                    val = genes
                elif metric == 'bases':
                    val = bases
                elif metric == 'genomes':
                    val = genomes
                else:
                    assert False
                metric2table[metric][KO].append(val)

    for metric in METRICS:
        metric2table[metric] = pd.DataFrame(metric2table[metric])

    # Calculate the percentage metrics
    for metric in METRICS:
        if metric == 'genes':
            s2norm = sdb.set_index('sample')['detected_genes'].to_dict()
        elif metric == 'bases':
            s2norm = sdb.set_index('sample')['bases_mapped_to_genes'].to_dict()
        elif metric == 'genomes':
            s2norm = sdb.set_index('sample')['detected_genomes'].to_dict()
        else:
            assert False

        db = metric2table[metric].copy()
        for d in TOTAL_KOS:
            db[d] = [x / s2norm[s] if s2norm[s] != 0 else 0 for x, s in zip(db[d], db['sample'])]

        metric2table[metric + '_fraction'] = db

    # Calculate the long-form table
    table = defaultdict(list)
    for sample, a2g2vals in s2a2g2vals.items():
        for a, g2vals in a2g2vals.items():
            for g, vals in g2vals.items():
                table['sample'].append(sample)
                table['anno'].append(a)
                table['genome'].append(g)
                for name, thing in zip(['genomes', 'genes', 'bases'], vals):
                    if name == 'genomes':
                        continue
                    table[name].append(thing)
    metric2table["long_data"] = pd.DataFrame(table)

    return metric2table

def load_s2a2g2vals(loc):
    """
    Create s2a2g2vals from long_data
    """

    s2a2g2vals = {}
    bdb = pd.read_csv(loc)
    for i, row in bdb.iterrows():
        s = row['sample']
        g = row['genome']
        a = row['anno']

        if s not in s2a2g2vals:
            s2a2g2vals[s] = {}
        if a not in s2a2g2vals[s]:
            s2a2g2vals[s][a] = {}
        s2a2g2vals[s][a][g] = [set([g]), row['genes'], row['bases']]

    return s2a2g2vals