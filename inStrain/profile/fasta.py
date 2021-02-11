'''
Opperations related to fasta files
'''

import logging
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

import inStrain.genomeUtilities

def load_fasta(**kwargs):
    '''
    Load the sequences to be profiled. All arguments from the "args" variable

    Return a table listing scaffold name, start, end
    '''
    fasta_loc = kwargs.get('fasta')
    window_length = kwargs.get('window_length')

    scaffolds_to_profile = kwargs.get('scaffolds_to_profile', None)
    use_full_fasta_header = kwargs.get('use_full_fasta_header', False)

    if use_full_fasta_header:
        scaff2sequence = {r.description:r.seq.upper() for r in SeqIO.parse(fasta_loc, "fasta")}
    else:
        scaff2sequence = {r.id:r.seq.upper() for r in SeqIO.parse(fasta_loc, "fasta")}

    # Get scaffold2length
    s2l = {s:len(scaff2sequence[s]) for s in list(scaff2sequence.keys())}

    # Generate splits
    table = defaultdict(list)
    for scaffold, sLen in s2l.items():
        for i, (split_start, split_end) in enumerate(iterate_splits(sLen, window_length)):
            table['scaffold'].append(scaffold)
            table['split_number'].append(i)
            table['start'].append(split_start)
            table['end'].append(split_end)
    Fdb = pd.DataFrame(table)

    # This just takes too long, sadly
    # _validate_splits(Fdb, s2l)

    if scaffolds_to_profile != None:
        Fdb = Fdb[Fdb['scaffold'].isin(scaffolds_to_profile)]
        s2s = {scaff:seq for scaff, seq in scaff2sequence.items() if scaff in scaffolds_to_profile}

    if len(Fdb) == 0:
        logging.error("The provided scaffold list has no overlap with the provided .fasta file!")
        logging.error("Example scaffolds in list: {0}".format("\n".join(scaffolds_to_profile)))
        sys.exit()

    return Fdb, scaff2sequence, s2l # also return s2l - alexcc 5/9/2019: Nah, make it scaff2sequence (s2s) (M.O. 6/10/19)

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

def filter_fasta(FAdb, s2p, s2l, rl, **kwargs):
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
    cdb.loc[:, 'read_pairs'] = cdb['scaffold'].map(s2p)
    cdb.loc[:, 'length'] = cdb['scaffold'].map(s2l)

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
