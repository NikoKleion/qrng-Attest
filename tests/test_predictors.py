# the array forms of the SP 800-90B predictors on bits against the per-sample loops translated from NIST's reference
import contextlib

import numpy as np

from qrng_attest import estimators as E
from qrng_attest import predictors as P


def _sources(n, seed):
    rng = np.random.default_rng([seed, n])
    yield "uniform", rng.integers(0, 2, n)
    yield "biased", (rng.random(n) < 0.9).astype(int)
    yield "zeros", np.zeros(n, int)
    yield "period 7", (np.arange(n) % 7 < 3).astype(int)
    yield "sticky", np.cumsum(rng.random(n) < 0.05) % 2
    s = rng.integers(0, 2, n)
    if n > 37:
        s[37:] = np.where(rng.random(n - 37) < 0.8, s[:-37], s[37:])
    yield "lag 37", s


@contextlib.contextmanager
def _set(**values):
    saved = {name: getattr(P, name) for name in values}
    try:
        for name, value in values.items():
            setattr(P, name, value)
        yield
    finally:
        for name, value in saved.items():
            setattr(P, name, value)


def _equal(lengths, seed):
    for n in lengths:
        for name, S in _sources(n, seed):
            for pred in P.COUNTS:
                got, want = P.COUNTS[pred](S, 2), P.LOOPS[pred](S, 2)
                assert got == want, (n, name, pred, got, want)


def test_array_forms_equal_the_loops_on_bits():
    # lengths on both sides of each threshold: 3 samples, the 16 symbols of LZ78Y, the 63 and 4095 sample windows
    _equal([3, 4, 17, 18, 19, 63, 64, 65, 129, 700], 1)
    _equal([4095, 4096, 4097, 5000], 2)


def test_array_forms_equal_the_loops_across_block_edges():
    for block, lengths in ((16, [18, 100, 400]), (100, [18, 100, 700]), (2048, [700, 4200])):
        with _set(BLOCK=block):
            _equal(lengths, 3)


def test_array_forms_equal_the_loops_where_the_limits_bind():
    # depth 6 holds 128 pairs against a limit of 100, and lengths 1 to 5 hold 62 contexts against a limit of 40
    with _set(D_MMC=6, MMC_MAX_ENTRIES=100, B_LZ=5, LZ_MAX_DICTIONARY=40):
        S = np.random.default_rng(4).integers(0, 2, 3000)
        ctx = sum(S[5 - d:len(S) - 1 - d] << d for d in range(6))
        assert len(np.unique(2 * ctx + S[6:])) > 100
        for block in (64, 1 << 19):
            with _set(BLOCK=block):
                _equal([18, 60, 300, 3000], 4)
    with _set(D_MMC=8, MMC_MAX_ENTRIES=300, B_LZ=7, LZ_MAX_DICTIONARY=9):
        _equal([60, 2000], 5)


def test_multimmc_takes_the_loop_where_the_limit_binds_above_the_last_depth():
    # 64 pairs at depth 4 against a limit of 40: outside what the array form covers, and it differs there
    with _set(D_MMC=6, MMC_MAX_ENTRIES=40):
        S = np.random.default_rng(6).integers(0, 2, 3000)
        assert P.multimmc_counts(S, 2) == P.multimmc_loop(S, 2)
        assert P._multimmc_bits(S) != P.multimmc_loop(S, 2)


def test_lag_array_form_on_wider_alphabets():
    rng = np.random.default_rng(7)
    for k in (3, 4, 16, 256, 1000):
        for n in (3, 50, 129, 2000):
            S = rng.integers(0, k, n)
            assert P.lag_counts(S) == P.lag_loop(S), (k, n)
            S = (S // 2) * 2 % k
            assert P.lag_counts(S) == P.lag_loop(S), (k, n)


def test_wider_alphabets_take_the_loops():
    S = np.random.default_rng(8).integers(0, 4, 5000)
    for pred in ("multimcw", "multimmc", "lz78y"):
        assert P.COUNTS[pred](S, 4) == P.LOOPS[pred](S, 4)


def test_longest_run_and_leaders():
    assert P._longest_run(np.array([], bool)) == 0 and P._longest_run(np.zeros(5, bool)) == 0
    assert P._longest_run(np.array([1, 1, 0, 1, 1, 1, 0, 1], bool)) == 3 and P._longest_run(np.ones(4, bool)) == 4
    # three predictors, five steps: a tie goes to the later row, and a hit below the best score leaves the leader
    hits = np.array([[1, 0, 1, 0, 0],
                     [1, 0, 0, 0, 1],
                     [0, 1, 0, 0, 1]], bool)
    prior, scores, leader = P._leaders(hits, np.zeros(3, np.int32), 0)
    assert prior.tolist() == [0, 1, 2, 0, 0] and scores.tolist() == [2, 2, 2] and leader == 2


def test_assessment_of_a_million_bits_holds_its_values():
    # every estimator on a million uniform bits, the size the standard asks for
    S = np.random.default_rng(0).integers(0, 2, 1000000)
    est = E.all_estimators(S, 2)
    assert len(est) == 10 and all(0.8 < v <= 1.0 for v in est.values()), est
    assert P.lag_counts(S) == (500112, 999999, 18) and P.multimcw_counts(S, 2) == (499919, 999937, 20)
    assert P.multimmc_counts(S, 2) == (501062, 999998, 18) and P.lz78y_counts(S, 2) == (499314, 999983, 18)
