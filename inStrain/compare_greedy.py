#!/usr/bin/env python

"""
None of this code is in use at the moment
"""

def greedy_main(RCprof, names, Sprofiles, scaffolds_to_compare, s2l, **kwargs):
    '''
    DEPRECATED

    Perform greedy clustering instead of all-vs-all comparisons
    '''
    g_ani = kwargs.get('g_ani', 0.99)
    g_cov = kwargs.get('g_cov', 0.5)
    g_mm = kwargs.get('g_mm', 100)

    # Set up the clustering
    name2cluster = {}
    cluster2rep = {}
    name2Sprofile = {n:s for n, s in zip(names, Sprofiles)}

    # Set up the initial cluster
    clusters = 0
    name2cluster[names[0]] = clusters
    cluster2rep[clusters] = names[0]
    clusters += 1

    # Make sure the names are unique
    assert len(names) == len(set(names)), "IS objects do not have unique names!"

    # Figure out the total bin length
    if kwargs.get('scaffolds', None) != None:
        scaffolds = inStrain.profile.fasta.load_scaff_list(kwargs.get('scaffolds', None))
        BIN_LENGTH = sum([s2l[s] for s in scaffolds])
    else:
        BIN_LENGTH = sum([l for s, l in s2l.items()])
    logging.info("BIN_LENGTH is {0}".format(BIN_LENGTH))

    # Do the iterations
    ndbs = []
    ddbs = []
    for name, IS in zip(names[1:], Sprofiles[1:]):

        # Get a list of Sprofiles to compare to
        To_compare = []
        compare_names = []
        for i in range(0, clusters):
            To_compare.append(name2Sprofile[cluster2rep[i]])
            compare_names.append(cluster2rep[i])

        # Get the distance matrix
        Ddb, Ndb = compare_Sprofiles_wrapper(IS, To_compare, name, compare_names, scaffolds_to_compare, s2l, BIN_LENGTH, **kwargs)
        ndbs.append(Ndb)
        ddbs.append(Ddb)

        # Adjust the clusters
        found = False
        for i, row in Ddb.iterrows():
            if (row['popANI'] >= g_ani) & (row['cov'] >= g_cov):
                found = True
                logging.debug("{0} is in the same cluster as {1}".format(name, row['name2']))
                name2cluster[name] = name2cluster[row['name2']]
                break

        if not found:
            logging.debug("{0} is a new cluster".format(name))
            name2cluster[name] = clusters
            cluster2rep[clusters] = name
            clusters += 1

    # Make the output
    Cdb = pd.DataFrame(list(name2cluster.items()), columns=['name', 'cluster'])
    Ndb = pd.concat(ndbs)
    Ddb = pd.concat(ddbs)

    # Store results
    outbase = RCprof.get_location('output') + os.path.basename(RCprof.get('location')) + '_'

    RCprof.store('greedy_clusters', Cdb, 'pandas', 'Cluster affiliations of the requested IS objects')
    RCprof.store('scaffold2length', s2l, 'dictionary', 'Scaffold to length')
    RCprof.store('comparisonsTable_greedy', Ndb, 'pandas', 'Comparisons between the requested IS objects done from a greedy clustering')
    RCprof.store('parsed_comparisonsTable_greedy', Ddb, 'pandas', 'Parsed comparisons between the requested IS objects done from a greedy clustering')

    Cdb.to_csv(outbase + 'greedyClusters.tsv', index=False, sep='\t')
    Ddb.to_csv(outbase + 'parsed_comparisonsTable_greedy.tsv', index=False, sep='\t')

def compare_Sprofiles_wrapper(IS1, IS_list, name1, names, scaffolds_to_compare,
                                s2l, BIN_LENGTH, **kwargs):
    '''
    Compare IS1 to every IS in the IS_list
    '''
    table = defaultdict(list)
    cdbs = []
    for cIS, name2 in zip(IS_list, names):
        results = compare_Sprofiles(IS1, cIS, [name1, name2], scaffolds_to_compare, s2l, BIN_LENGTH, **kwargs)
        Cdb, ANI, cov = results
        table['name1'].append(name1)
        table['name2'].append(name2)
        table['ANI'].append(ANI)
        table['cov'].append(cov)
        for thing in ['g_ani', 'g_cov', 'g_mm']:
            table[thing].append(kwargs.get(thing, np.nan))
        cdbs.append(Cdb)

    Ndb = pd.concat(cdbs)
    Ddb = pd.DataFrame(table)
    return Ddb, Ndb

def compare_Sprofiles(IS1, cIS, names, scaffolds_to_compare, s2l, BIN_LENGTH, **kwargs):
    '''
    Compare all scaffolds of two Sprofiles
    '''
    Cdb, scaff2pair2mm2SNPs, scaff2pair2mm2cov = compare_scaffolds(names, [IS1, cIS], scaffolds_to_compare, s2l, **kwargs)
    ANI, cov = calc_cov_ani(Cdb, BIN_LENGTH, **kwargs)
    return Cdb, ANI, cov

def calc_cov_ani(Cdb, BIN_LENGTH, **kwargs):
    '''
    Return ANI and coverage assuming all scaffolds belong to the same genome
    '''
    g_ani = kwargs.get('g_ani', 0.99)
    g_cov = kwargs.get('g_cov', 0.5)
    g_mm = kwargs.get('g_mm', 100)

    db = Cdb[Cdb['mm'] <= g_mm].sort_values('mm').drop_duplicates(subset=['scaffold', 'name1', 'name2'], keep='last')
    tcb = sum(db['compared_bases_count'])

    if tcb == 0:
        popANI = np.nan
    else:
        popANI = sum(a * c  if a == a else 0 for a, c in
                            zip(db['popANI'], db['compared_bases_count'])) / tcb

    return popANI, (tcb / BIN_LENGTH)