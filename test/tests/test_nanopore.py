"""
Tests for running inStrain on error-prone long reads (e.g. Nanopore)

These simulate long, single-end reads with ~3% substitution errors and some indels from a mix of two
strains, so the expected SNVs are known exactly. They check that:

* with the default settings (Q30, which Nanopore bases rarely reach) almost nothing is counted
* lowering --min_base_quality without a matching error model calls many false SNVs
* lowering --min_base_quality with the matching error model calls the real SNVs and (almost) nothing else

Long reads are profiled with --skip_mm_profiling: each read has dozens of mismatches, so mm-level profiling is slow
and makes real SNVs flicker in and out at low mm levels, which marks them as 'cryptic'.
"""

import importlib
import logging
import os
import random

import pysam
import pytest

import inStrain
import inStrain.argumentParser
import inStrain.controller
import inStrain.profile.snv_utilities
import inStrain.SNVprofile

BASES = 'ACGT'


def _mutate(base, rng):
    return rng.choice([b for b in BASES if b != base])


def simulate_long_reads(out_dir, seed=42, genome_length=20000, num_reads=600, read_length=(1500, 3000),
                        sub_rate=0.03, ins_rate=0.005, del_rate=0.005, base_quality=15,
                        num_snvs=25, strain_fraction=0.5):
    """
    Simulate a reference, a second strain carrying num_snvs SNVs, and error-prone single-end reads from a
    mix of the two strains. Reads are written directly as alignments (no aligner needed).

    Returns (fasta, bam, snv_positions) where snv_positions are the 0-based true SNV positions
    """
    rng = random.Random(seed)
    scaffold = 'sim_scaffold'
    reference = ''.join(rng.choice(BASES) for _ in range(genome_length))

    # Second strain
    snv_positions = sorted(rng.sample(range(500, genome_length - 500), num_snvs))
    strain2 = list(reference)
    for p in snv_positions:
        strain2[p] = _mutate(reference[p], rng)
    strain2 = ''.join(strain2)

    fasta = os.path.join(out_dir, 'sim.fasta')
    with open(fasta, 'w') as o:
        o.write(f'>{scaffold}\n{reference}\n')

    header = {'HD': {'VN': '1.6', 'SO': 'unsorted'}, 'SQ': [{'SN': scaffold, 'LN': genome_length}]}
    unsorted_bam = os.path.join(out_dir, 'sim.unsorted.bam')
    with pysam.AlignmentFile(unsorted_bam, 'wb', header=header) as out:
        for r in range(num_reads):
            source = strain2 if rng.random() < strain_fraction else reference
            length = rng.randint(*read_length)
            start = rng.randint(0, genome_length - length)

            seq, cigar, nm = [], [], 0
            ref_pos = start

            def add(op, n=1):
                if cigar and cigar[-1][0] == op:
                    cigar[-1] = (op, cigar[-1][1] + n)
                else:
                    cigar.append((op, n))

            while ref_pos < start + length:
                x = rng.random()
                if x < del_rate and seq:
                    add(2)  # deletion
                    nm += 1
                    ref_pos += 1
                    continue
                if x < del_rate + ins_rate and seq:
                    seq.append(rng.choice(BASES))
                    add(1)  # insertion
                    nm += 1
                    continue

                base = source[ref_pos]
                if rng.random() < sub_rate:
                    base = _mutate(base, rng)
                if base != reference[ref_pos]:
                    nm += 1
                seq.append(base)
                add(0)  # match / mismatch
                ref_pos += 1

            a = pysam.AlignedSegment()
            a.query_name = f'read_{r}'
            a.query_sequence = ''.join(seq)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = cigar
            a.query_qualities = pysam.qualitystring_to_array(chr(base_quality + 33) * len(seq))
            a.set_tag('NM', nm)
            out.write(a)

    bam = os.path.join(out_dir, 'sim.bam')
    pysam.sort('-o', bam, unsorted_bam)
    pysam.index(bam)
    os.remove(unsorted_bam)

    return fasta, bam, snv_positions


@pytest.fixture(scope='module')
def sim(tmp_path_factory):
    out_dir = str(tmp_path_factory.mktemp('nanopore_sim'))
    fasta, bam, snv_positions = simulate_long_reads(out_dir)
    return {'dir': out_dir, 'fasta': fasta, 'bam': bam, 'snvs': set(snv_positions)}


def _profile(sim, name, extra_args):
    out = os.path.join(sim['dir'], name)
    importlib.reload(logging)  # so inStrain sets up its own log file (pytest installs a log handler)
    cmd = (f"profile {sim['bam']} {sim['fasta']} -o {out} -p 1 --skip_plot_generation "
           f"--pairing_filter non_discordant -l 0.9 {extra_args}")
    inStrain.controller.Controller().main(inStrain.argumentParser.parse_args(cmd.split()))
    IS = inStrain.SNVprofile.SNVprofile(out)

    sdb = IS.get_nonredundant_snv_table()
    called = set() if len(sdb) == 0 else set(sdb[sdb['allele_count'] >= 2]['position'].astype(int))
    scaffold_table = IS.get_nonredundant_scaffold_table()
    coverage = 0 if len(scaffold_table) == 0 else scaffold_table['coverage'].iloc[0]
    return called, coverage


def test_nanopore_default_settings_count_nothing(sim):
    """
    With the default Q30 cutoff the Q15 bases are not counted at all (the issue reported in #99)
    """
    called, coverage = _profile(sim, 'default', '')
    assert coverage < 1, coverage
    assert len(called) == 0, called

    # The user should be told how to fix it
    with open(os.path.join(sim['dir'], 'default', 'log', 'log.log')) as o:
        log = o.read()
    assert 'looks like long reads' in log
    assert '--min_base_quality' in log and '--skip_mm_profiling' in log


def test_nanopore_needs_matching_error_model(sim):
    """
    Counting Q15 bases but keeping the Illumina (Q30) error model calls many false SNVs
    """
    called, coverage = _profile(sim, 'wrong_model', '--min_base_quality 15 --error_rate 0.001 --skip_mm_profiling')
    assert coverage > 50, coverage
    false_positives = called - sim['snvs']
    assert len(false_positives) > 20, len(false_positives)


def test_nanopore_with_matching_error_model(sim):
    """
    Counting Q15 bases with the matching error model (calculated from --min_base_quality) finds the real
    SNVs and essentially nothing else
    """
    called, coverage = _profile(sim, 'right_model', '--min_base_quality 15 --skip_mm_profiling')
    assert coverage > 50, coverage

    found = called & sim['snvs']
    false_positives = called - sim['snvs']
    assert len(found) >= len(sim['snvs']) - 1, (len(found), len(sim['snvs']))
    assert len(false_positives) <= 2, sorted(false_positives)


def test_error_rate_model_matches_shipped_model():
    """
    The calculated model at a Q30 error rate should agree with the simulated table shipped with inStrain.
    The table is a Monte Carlo simulation, so allow a difference of one read up to 1000x coverage and a few
    reads beyond that (where the simulation runs out of resolution)
    """
    table = inStrain.profile.snv_utilities.load_null_model()
    calculated = inStrain.profile.snv_utilities.generate_error_rate_model(0.001)
    shared = [c for c in table if c > 0 and c in calculated]
    assert len(shared) > 9000
    assert all(abs(table[c] - calculated[c]) <= 1 for c in shared if c <= 1000)
    assert all(abs(table[c] - calculated[c]) <= 4 for c in shared)


def test_load_null_model_options():
    default = inStrain.profile.snv_utilities.load_null_model()
    legacy = inStrain.profile.snv_utilities.load_null_model(legacy_snv_thresholds=True)
    q20 = inStrain.profile.snv_utilities.load_null_model(min_base_quality=20)
    q20_explicit = inStrain.profile.snv_utilities.load_null_model(error_rate=0.01)

    assert all(legacy[c] == default[c] - 1 for c in default if c > 0)
    assert q20 == q20_explicit
    # Higher error rates need more supporting reads
    assert all(q20[c] >= default[c] for c in [10, 50, 100, 1000])

    with pytest.raises(ValueError):
        inStrain.profile.snv_utilities.generate_error_rate_model(0.9)
