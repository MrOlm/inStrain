#!/usr/bin/env python3
"""
Re-apply stricter SNV calling thresholds (--min_cov / --min_freq) to an existing inStrain profile,
without re-running inStrain profile.

Raising the thresholds can only remove calls, and an inStrain profile already stores everything needed
to recalculate them:

* raw_snp_table: every site called at the original thresholds, with its A/C/T/G counts at each
  read-mismatch (mm) level. These sites are re-called with the new thresholds.
* covT: the coverage of every position at each mm level. Used for breadth_minCov and to find positions
  that are no longer counted.
* clonT: the nucleotide diversity of every position with at least the original min_cov. Nucleotide
  diversity doesn't depend on min_freq, so positions below the new min_cov are just removed.

The scaffold- and genome-level tables are then rebuilt with inStrain's own code, so the results match
running inStrain profile with the new thresholds. The exceptions are linkage (which needs the reads, so it's
left out) and gene-level results (re-run them with -g, as for inStrain profile_genes).

The original profile is not changed; the re-filtered profile is written to a new folder.

Example:
    python refilter_profile.py -i sample.IS -o sample.c10_f10.IS --min_cov 10 --min_freq 0.1
"""

__author__ = "Matt Olm"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import logging
import os
import shutil
import sys
from collections import defaultdict
import re

import numpy as np
import pandas as pd

import inStrain
import inStrain.controller
import inStrain.genomeUtilities
import inStrain.GeneProfile
import inStrain.profile.profile_utilities as profile_utilities
import inStrain.profile.snv_utilities as snv_utilities
import inStrain.SNVprofile


SNV_COLUMNS = ['scaffold', 'position', 'ref_base', 'A', 'C', 'T', 'G', 'con_base', 'var_base', 'mm',
               'allele_count', 'class', 'cryptic', 'position_coverage']


def recall_scaffold(raw, covT, null_model, min_cov, min_freq):
    """
    Re-call the SNVs of one scaffold with new thresholds

    This follows inStrain.profile.snv_utilities.update_snp_table: each position is evaluated at each mm level
    with reads, and a position is "cryptic" if it was an SNV at a lower mm level but not at a higher one.

    Args:
        raw: raw_snp_table rows for this scaffold (sites called at the original, lower, thresholds)
        covT: covT for this scaffold (mm -> per-position read counts at exactly that mm)

    Returns:
        DataFrame in the format of raw_snp_table
    """
    if (raw is None) or (len(raw) == 0):
        return pd.DataFrame(columns=SNV_COLUMNS)

    positions = set(raw['position'].astype(int))
    pos2ref = dict(zip(raw['position'].astype(int), raw['ref_base']))
    rows = {(int(p), int(mm)): [a, c, t, g] for p, mm, a, c, t, g in
            zip(raw['position'], raw['mm'], raw['A'], raw['C'], raw['T'], raw['G'])}

    # For each position, the read count at each mm level that has reads
    pos2mm2count = defaultdict(dict)
    for mm, counts in covT.items():
        counts = counts[counts.index.isin(positions)]
        for p, c in zip(counts.index, counts.values):
            if c > 0:
                pos2mm2count[int(p)][int(mm)] = c

    table = defaultdict(list)
    p2c = {}
    for p in sorted(positions):
        ref_base = pos2ref[p]
        anySNP = False
        coverage = 0
        for mm in sorted(pos2mm2count[p].keys()):
            coverage += pos2mm2count[p][mm]

            if (p, mm) in rows:
                counts = rows[(p, mm)]
                snp, morphia = snv_utilities.call_snv_site(counts, ref_base, null_model,
                                                           min_cov=min_cov, min_freq=min_freq)
            elif coverage >= min_cov:
                # Not called originally despite enough coverage, so only the reference base was present
                snp, morphia = -1, 1
            else:
                snp, morphia = None, 0

            if not snp:
                continue

            elif snp != -1:
                counts_temp = list(counts)
                counts_temp[snv_utilities.P2C[snp]] = 0
                var_base = snv_utilities.C2P[counts_temp.index(sorted(counts_temp)[-1])]
                snp_class = snv_utilities.calc_snp_class(snp, ref_base, var_base, counts, morphia, null_model,
                                                         min_cov=min_cov, min_freq=min_freq)

                table['position'].append(p)
                table['ref_base'].append(ref_base)
                for b, c in zip(['A', 'C', 'T', 'G'], counts):
                    table[b].append(c)
                table['con_base'].append(snp)
                table['var_base'].append(var_base)
                table['mm'].append(mm)
                table['allele_count'].append(morphia)
                table['class'].append(snp_class)

                if morphia >= 2:
                    anySNP = True
                elif (morphia == 1) & anySNP:
                    p2c[p] = True

            elif anySNP:
                p2c[p] = True

    db = pd.DataFrame(table)
    if len(db) == 0:
        return pd.DataFrame(columns=SNV_COLUMNS)
    db['scaffold'] = raw['scaffold'].iloc[0]
    db['cryptic'] = db['position'].map(p2c).fillna(False).astype(bool)
    db['position_coverage'] = db['A'] + db['C'] + db['T'] + db['G']
    return db[SNV_COLUMNS]


def mask_clonT(clonT, covT, min_cov):
    """
    Remove positions with less than min_cov coverage from clonT (scaffold-level)
    """
    masked = {}
    for mm, clons in clonT.items():
        coverage = profile_utilities.mm_counts_to_counts_shrunk(covT, mm)
        coverage = coverage.reindex(clons.index).fillna(0)
        masked[mm] = clons[(coverage >= min_cov).values]
    return masked


def _version_tuple(version):
    """
    '1.3.0g' -> (1, 3, 0)
    """
    return tuple(int(x) for x in re.findall(r'\d+', str(version))[:3])


def load_null_model(IS, args):
    """
    Use the same null model that the profile was made with. Profiles from before v1.12 used thresholds that
    were one read lower (see --legacy_snv_thresholds in inStrain profile)
    """
    legacy = args.legacy_snv_thresholds
    if legacy is None:
        legacy = _version_tuple(IS.get('version')) < (1, 12)
        if legacy:
            logging.info(f"This profile was made with inStrain v{IS.get('version')}, so the SNV thresholds of "
                         "inStrain <1.12 are used to match it (override with --current_snv_thresholds)")

    return snv_utilities.load_null_model(fdr=args.fdr, min_base_quality=args.min_base_quality,
                                         error_rate=args.error_rate, legacy_snv_thresholds=legacy)


def refilter_profile(args):
    if (args.min_cov < args.original_min_cov) or (args.min_freq < args.original_min_freq):
        raise ValueError("The new --min_cov / --min_freq can't be lower than the ones the profile was made with "
                         "(--original_min_cov / --original_min_freq); the sites that would be gained are not stored")
    if os.path.exists(args.output):
        raise ValueError(f"{args.output} already exists")

    # Copy the profile and set up logging
    shutil.copytree(args.input, args.output)
    IS = inStrain.SNVprofile.SNVprofile(args.output)
    for f in os.listdir(IS.get_location('log')):
        os.remove(os.path.join(IS.get_location('log'), f))
    inStrain.controller.setup_logger(IS.get_location('log') + 'log.log')
    logging.info(f"Re-filtering {args.input} with --min_cov {args.min_cov} --min_freq {args.min_freq}")

    # Output files are re-made below (and named after the new profile)
    for f in os.listdir(IS.get_location('output')):
        os.remove(os.path.join(IS.get_location('output'), f))

    null_model = load_null_model(IS, args)

    raw = IS.get('raw_snp_table')
    covT = IS.get('covT')
    clonT = IS.get('clonT')
    clonTR = IS.get('clonTR')
    s2l = IS.get('scaffold2length')
    scaffolds = list(IS.get('cumulative_scaffold_table')['scaffold'].astype(str).unique())

    scaff2raw = {}
    if (raw is not None) and (len(raw) > 0):
        raw = raw.copy()
        raw['scaffold'] = raw['scaffold'].astype(str)
        scaff2raw = {s: db for s, db in raw.groupby('scaffold', observed=True)}

    new_raw, new_scaffold_tables, new_clonT = [], [], {}
    for scaffold in scaffolds:
        sdb = recall_scaffold(scaff2raw.get(scaffold), covT[scaffold], null_model, args.min_cov, args.min_freq)
        new_raw.append(sdb)
        new_clonT[scaffold] = mask_clonT(clonT[scaffold], covT[scaffold], args.min_cov)
        new_scaffold_tables.append(profile_utilities.make_coverage_table(
            covT[scaffold], new_clonT[scaffold], clonTR[scaffold], s2l[scaffold], scaffold,
            sdb, null_model, min_freq=args.min_freq))

    new_raw = pd.concat(new_raw).reset_index(drop=True)
    for col in ['A', 'C', 'G', 'T', 'mm', 'position']:
        new_raw[col] = new_raw[col].astype(int)
    cumulative_snv_table = profile_utilities._parse_Sdb(profile_utilities._make_snp_table(new_raw.copy()))

    IS.store('raw_snp_table', new_raw, 'pandas', 'Contains raw SNP information on a mm level')
    IS.store('cumulative_snv_table', cumulative_snv_table, 'pandas', 'Cumulative SNP on mm level')
    IS.store('cumulative_scaffold_table', pd.concat(new_scaffold_tables).reset_index(drop=True), 'pandas',
             'Cumulative coverage on mm level')
    IS.store('clonT', new_clonT, 'special', 'Scaffold -> mm -> position based clonality')

    # Linkage needs the reads, so it can't be re-filtered
    IS.store('raw_linkage_table', pd.DataFrame(), 'pandas', 'Raw table of linkage information')
    logging.info("Linkage is not re-calculated (it needs the reads), so it is left out of the output")

    # Genes
    if args.gene_file is not None:
        inStrain.GeneProfile.Controller().main(argparse.Namespace(
            IS=args.output, gene_file=args.gene_file, processes=args.processes))
    elif IS.get('genes_table') is not None:
        logging.warning("This profile has gene-level results, which are not updated; pass the genes file "
                        "with -g to re-calculate them")

    # Genome-level information, using the scaffold-to-bin information stored in the profile
    if IS.get('scaffold2bin') is not None:
        GIdb = inStrain.genomeUtilities.genomeLevel_from_IS(IS)
        IS.store('genome_level_info', GIdb, 'pandas', 'Table of genome-level information')

    for t in ['SNVs', 'scaffold_info', 'genome_info', 'mapping_info'] + \
             (['gene_info'] if args.gene_file is not None else []):
        IS.generate(t)

    logging.info(f"Done! Re-filtered profile written to {args.output}")


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='inStrain profile to re-filter')
    parser.add_argument('-o', '--output', required=True, help='Where to write the re-filtered profile (must not exist)')
    parser.add_argument('-c', '--min_cov', type=int, required=True, help='New minimum coverage')
    parser.add_argument('-f', '--min_freq', type=float, required=True, help='New minimum SNV frequency')
    parser.add_argument('-g', '--gene_file', default=None,
                        help='Genes file, to re-calculate gene-level results (as for inStrain profile_genes)')
    parser.add_argument('-p', '--processes', type=int, default=6, help='Processes to use for gene profiling')

    o = parser.add_argument_group('Settings the profile was made with (defaults are the inStrain profile defaults)')
    o.add_argument('--original_min_cov', type=int, default=5)
    o.add_argument('--original_min_freq', type=float, default=0.05)
    o.add_argument('--fdr', type=float, default=1e-6)
    o.add_argument('--min_base_quality', type=int, default=30)
    o.add_argument('--error_rate', type=float, default=None)
    o.add_argument('--legacy_snv_thresholds', dest='legacy_snv_thresholds', action='store_true', default=None,
                   help='Use the SNV thresholds of inStrain <1.12 (default: if the profile was made with inStrain <1.12)')
    o.add_argument('--current_snv_thresholds', dest='legacy_snv_thresholds', action='store_false',
                   help='Use the current SNV thresholds, even for profiles made with inStrain <1.12')
    return parser.parse_args(args)


if __name__ == '__main__':
    refilter_profile(parse_args())
