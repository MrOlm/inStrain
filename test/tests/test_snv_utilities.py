"""
Tests for SNV utility functions
"""

import os
import tempfile

import inStrain
import inStrain.profile.snv_utilities


def _null_model_loc():
    """The NullModel.txt that ships with inStrain, resolved as the controllers do."""
    return os.path.join(os.path.dirname(inStrain.__file__), 'helper_files', 'NullModel.txt')


def test_generate_snp_model_stores_read_counts():
    """The model must store the read count whose probability is below fdr, not its index.

    Uses a synthetic table so the assertion depends only on the parsing logic, not on
    the Monte-Carlo values in the shipped file. With fdr=1e-6:
      coverage 10 -> the first column below 1e-6 is the 3rd (k=3)
      coverage 20 -> the first column below 1e-6 is the 2nd (k=2)
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        table = os.path.join(tmpdir, 'SyntheticNullModel.txt')
        with open(table, 'w') as o:
            o.write('coverage\t1\t2\t3\t4\n')
            o.write('10\t0.01\t0.001\t0.0000001\t0.0\n')
            o.write('20\t0.02\t0.0000005\t0.0\t0.0\n')

        model = inStrain.profile.snv_utilities.generate_snp_model(table, fdr=1e-6)
        assert model[10] == 3, model[10]
        assert model[20] == 2, model[20]
        assert model[-1] == 3, model[-1]

        legacy = inStrain.profile.snv_utilities.generate_snp_model(
            table, fdr=1e-6, legacy_thresholds=True)
        assert legacy[10] == 2, legacy[10]
        assert legacy[20] == 1, legacy[20]


def test_generate_snp_model_on_shipped_table():
    """Spot-check the real table.

    At coverage 30 the row is 0.02998 (k=1), 0.00018 (k=2), 1e-06 (k=3), 0.0 (k=4);
    the first value strictly below 1e-6 is k=4, so an allele needs 4 reads.
    At coverage 5 the row reaches 0.0 at k=3.
    """
    model = inStrain.profile.snv_utilities.generate_snp_model(_null_model_loc(), fdr=1e-6)
    assert model[5] == 3, model[5]
    assert model[30] == 4, model[30]

    legacy = inStrain.profile.snv_utilities.generate_snp_model(
        _null_model_loc(), fdr=1e-6, legacy_thresholds=True)
    assert legacy[5] == 2, legacy[5]
    assert legacy[30] == 3, legacy[30]

    # legacy_thresholds must reproduce the previous behaviour at every coverage.
    assert all(legacy[c] == model[c] - 1 for c in model if c > 0)
