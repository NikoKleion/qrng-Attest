# SP 800-90B min-entropy estimators (section 6.3) and the multi-bit combination
import math
import numpy as np

from . import predictors as _P

Z = 2.5758293035489008


def _hmin(p):
    return float(-math.log2(min(max(p, 1e-12), 1.0)))


def _upper_bound_p(phat, n):
    # upper 99% confidence bound on a probability from n samples
    if n <= 1:
        return min(1.0, phat)
    return float(min(1.0, phat + Z * math.sqrt(phat * (1.0 - phat) / (n - 1))))


# 6.3.1 Most Common Value
def most_common_value(S, k=None):
    L = len(S)
    counts = np.bincount(S, minlength=(k or (int(S.max()) + 1)))
    phat = counts.max() / L
    return _hmin(_upper_bound_p(phat, L))


# 6.3.2 Collision (binary)
def collision(S, k=2):
    # collision estimate, binary only. Uses the mean number of samples until a value repeats.
    if k != 2:
        return None
    times = []
    i = 0; n = len(S)
    while i < n:
        seen = set(); j = i; repeat = False
        while j < n:
            v = int(S[j]); j += 1
            if v in seen:
                repeat = True; break
            seen.add(v)
        if repeat:
            times.append(j - i)
        i = j
    if len(times) < 2:
        return 1.0
    t = np.array(times, float)
    mean = t.mean(); sd = t.std(ddof=1); v = len(t)
    mean_lb = mean - Z * sd / math.sqrt(v)
    disc = 5.0 - 2.0 * mean_lb
    if disc <= 0:
        return 1.0
    p = 0.5 + 0.5 * math.sqrt(disc)
    p = min(max(p, 0.5), 1.0)
    return _hmin(p)


# 6.3.3 Markov
def markov(S, k=None, path_len=128):
    # Markov estimate: entropy of the most likely path. Change path_len to set the path length.
    k = k or (int(S.max()) + 1)
    L = len(S)
    init = np.bincount(S, minlength=k).astype(float); init /= init.sum()
    trans = np.full((k, k), 0.0)
    for a, b in zip(S[:-1], S[1:]):
        trans[a, b] += 1.0
    rows = trans.sum(axis=1, keepdims=True)
    trans = np.divide(trans, rows, out=np.full_like(trans, 1.0 / k), where=rows > 0)
    logp = np.log2(np.clip(init, 1e-300, None))
    lt = np.log2(np.clip(trans, 1e-300, None))
    for _ in range(path_len - 1):
        logp = (logp[:, None] + lt).max(axis=0)
    best = float(logp.max())
    return float(min(math.log2(k), -best / path_len))


# 6.3.5 / 6.3.6 t-Tuple and LRS
def _prefix_ranks(S):
    # ranks[m][i] orders the substrings of length 2^m by position i; the last level orders the suffixes
    n = len(S)
    rank = np.unique(S, return_inverse=True)[1].astype(np.int64)
    ranks = [rank.astype(np.int32)]
    step = 1
    while int(rank.max()) < n - 1 and step < n:
        nxt = np.zeros(n, np.int64)
        nxt[:n - step] = rank[step:] + 1
        rank = np.unique(rank * (n + 1) + nxt, return_inverse=True)[1].astype(np.int64)
        ranks.append(rank.astype(np.int32))
        step *= 2
    return ranks


def _common_prefixes(ranks, n):
    # length of the common prefix of each pair of suffixes adjacent in sorted order
    order = np.argsort(ranks[-1], kind="stable")
    a, b = order[:-1].astype(np.int64), order[1:].astype(np.int64)
    lcp = np.zeros(n - 1, np.int64)
    for m in range(len(ranks) - 1, -1, -1):
        inside = np.flatnonzero((a < n) & (b < n))
        same = inside[ranks[m][a[inside]] == ranks[m][b[inside]]]
        lcp[same] += 1 << m
        a[same] += 1 << m
        b[same] += 1 << m
    return lcp


_DENSE_LEVELS = 24


def _tuple_profile(S, k):
    # for each tuple length, the most common count and the collision count. Shared by t_tuple and lrs.
    S = np.asarray(S)
    n = len(S)
    if n < 2:
        return [(1, n, 0.0, n)] if n else []
    lcp = _common_prefixes(_prefix_ranks(S), n)
    v = int(lcp.max())
    most = np.ones(v + 2, np.int64)
    pairs = [0] * (v + 2)
    dense = min(v, _DENSE_LEVELS)
    for t in range(1, dense + 1):
        edge = np.diff(np.concatenate(([0], (lcp >= t).astype(np.int8), [0])))
        runs = np.flatnonzero(edge == -1) - np.flatnonzero(edge == 1)
        if runs.size:
            most[t] = int(runs.max()) + 1
            pairs[t] = int((runs * (runs + 1) // 2).sum())
    if v > dense:
        # above the dense levels, join neighbouring groups of suffixes from the longest common prefix down
        where = np.flatnonzero(lcp > dense)
        where = where[np.argsort(-lcp[where], kind="stable")]
        depth = lcp[where].tolist()
        where = where.tolist()
        first, last = {}, {}
        biggest, joined, i = 1, 0, 0
        for t in range(v, dense, -1):
            while i < len(where) and depth[i] >= t:
                j = where[i]
                lo, hi = first.get(j, j), last.get(j + 1, j + 1)
                joined += (j - lo + 1) * (hi - j)
                if hi - lo + 1 > biggest:
                    biggest = hi - lo + 1
                first[hi], last[lo] = lo, hi
                i += 1
            most[t], pairs[t] = biggest, joined
    prof = [(t, int(most[t]), float(pairs[t]), n - t + 1) for t in range(1, v + 1)]
    if v + 1 <= n:
        prof.append((v + 1, 1, 0.0, n - v))
    return prof


def t_tuple(S, k=None, _prof=None):
    L = len(S); k = k or int(S.max()) + 1
    prof = _prof if _prof is not None else _tuple_profile(S, k)
    best_p = 0.0
    for t, maxc, _coll, length in prof:
        if maxc < 35:
            break
        best_p = max(best_p, (maxc / length) ** (1.0 / t))
    if best_p <= 0:
        return 1.0
    return _hmin(_upper_bound_p(best_p, L))


def lrs(S, k=None, _prof=None):
    L = len(S); k = k or int(S.max()) + 1
    prof = _prof if _prof is not None else _tuple_profile(S, k)
    u = next((t for t, maxc, _c, _l in prof if maxc < 35), 1)
    w = max((t for t, maxc, _c, _l in prof if maxc >= 2), default=0)
    if w < u:
        return 1.0
    best_p = 0.0
    for t, _maxc, coll, length in prof:
        if u <= t <= w:
            denom = length * (length - 1) / 2.0
            if denom > 0 and coll > 0:
                best_p = max(best_p, (coll / denom) ** (1.0 / t))
    if best_p <= 0:
        return 1.0
    return _hmin(_upper_bound_p(best_p, L))


# 6.3.4 Compression
_COMP_B = 6
_COMP_ALPH = 1 << _COMP_B
_COMP_D = 1000


def _compression_G(z, d, num_blocks):
    # helper for the compression estimate: expected value for a symbol of probability z.
    v = num_blocks - d
    omz = 1.0 - z
    if omz <= 0.0:
        return 0.0
    i = np.arange(2, num_blocks + 1, dtype=np.float64)
    Bi = np.exp((i - 1.0) * math.log(omz))
    ai = np.log2(i) * Bi
    Ad1 = float(np.sum(ai[:d - 1]))
    Ai_minus_Ad1 = float(np.sum(ai[d - 1:]))
    tail_i = i[d - 1:num_blocks - 2]
    tail_ai = ai[d - 1:num_blocks - 2]
    first_sum = float(np.sum((num_blocks - tail_i) * tail_ai)) + (num_blocks - d) * Ad1
    return (1.0 / v) * z * (z * first_sum + Ai_minus_Ad1)


def _compression_expect(p, num_blocks):
    q = (1.0 - p) / (_COMP_ALPH - 1.0)
    return _compression_G(p, _COMP_D, num_blocks) + (_COMP_ALPH - 1.0) * _compression_G(q, _COMP_D, num_blocks)


def compression(S, k=2):
    # compression estimate, binary only. Returns bits per bit, or None if the stream is too short.
    if k != 2:
        return None
    bits = np.asarray(S, np.int64)
    n = len(bits)
    num_blocks = n // _COMP_B
    if num_blocks - _COMP_D < 2:
        return None
    blocks = np.zeros(num_blocks, np.int64)
    for j in range(_COMP_B):
        blocks = (blocks << 1) | bits[j:num_blocks * _COMP_B:_COMP_B]
    order = np.lexsort((np.arange(num_blocks), blocks))
    sb = blocks[order]
    prev = np.full(num_blocks, -1, np.int64)
    same = sb[1:] == sb[:-1]
    prev[order[1:][same]] = order[:-1][same]
    v = num_blocks - _COMP_D
    i = np.arange(_COMP_D, num_blocks)
    last1 = np.where(prev[_COMP_D:] >= 0, prev[_COMP_D:] + 1, 0)
    logs = np.log2((i + 1) - last1)
    x_bar = logs.mean()
    sigma = 0.5907 * math.sqrt(math.fsum(logs * logs) / (v - 1.0) - x_bar * x_bar)
    x_bar_prime = x_bar - Z * sigma / math.sqrt(v)
    if _compression_expect(1.0 / _COMP_ALPH, num_blocks) <= x_bar_prime:
        return 1.0
    lo, hi = 1.0 / _COMP_ALPH, 1.0
    p = 0.5 * (lo + hi)
    for _ in range(1076):
        val = _compression_expect(p, num_blocks)
        if abs(val - x_bar_prime) <= 4.0 * np.spacing(max(abs(val), abs(x_bar_prime))):
            break
        if x_bar_prime < val:
            lo = p
        else:
            hi = p
        new_p = 0.5 * (lo + hi)
        if new_p == p:
            break
        p = new_p
    if p <= 1.0 / _COMP_ALPH:
        return 1.0
    return float(-math.log2(p) / _COMP_B)


# 6.3.7-6.3.10 predictor framework
def _predictor_entropy(correct, k=2):
    # min-entropy from a record of which predictions were correct
    correct = np.asarray(correct, bool)
    n = len(correct)
    if n < 2:
        return math.log2(k)
    best = run = 0
    for hit in correct:
        run = run + 1 if hit else 0
        if run > best:
            best = run
    return _P.prediction_estimate(int(correct.sum()), n, best, k)


def multimcw(S, k=2):
    return _P.estimate("multimcw", S, k)


def lag(S, k=None):
    S = np.asarray(S)
    return _P.estimate("lag", S, k or int(S.max()) + 1)


def multimmc(S, k=2):
    return _P.estimate("multimmc", S, k)


def lz78y(S, k=2):
    return _P.estimate("lz78y", S, k)


def all_estimators(S, k=None, include_predictors=True):
    # run every estimator that applies to alphabet size k. Collision, Markov, and Compression are binary
    S = np.asarray(S).astype(int)
    if S.size == 0:
        return {}
    if S.min() < 0:
        raise ValueError("samples must be non-negative integers")
    k = max(int(k or 2), int(S.max()) + 1)
    prof = _tuple_profile(S, k)
    est = {"most_common_value": most_common_value(S, k),
           "t_tuple": t_tuple(S, k, _prof=prof),
           "lrs": lrs(S, k, _prof=prof)}
    if k == 2:
        est["collision"] = collision(S, k)
        est["markov"] = markov(S, k)
        est["compression"] = compression(S, k)
    if include_predictors:
        est.update(multimcw=multimcw(S, k), lag=lag(S, k), multimmc=multimmc(S, k), lz78y=lz78y(S, k))
    return {name: v for name, v in est.items() if v is not None}


def _to_bits(sym, word):
    # split each symbol into `word` bits, most significant first, and concatenate.
    sym = np.asarray(sym, np.int64)
    shifts = word - 1 - np.arange(word)
    return ((sym[:, None] >> shifts[None, :]) & 1).reshape(-1).astype(int)


def min_entropy(S, k=None, include_predictors=True, bits_per_symbol=None, relabel_bits=False):
    # main entry point. Binary sources are assessed directly. Multi-bit sources are assessed as symbols
    S = np.asarray(S).astype(int)
    if S.size and S.min() < 0:
        raise ValueError("samples must be non-negative integers")
    vals = np.unique(S)
    alph = int(len(vals))
    if alph <= 2:
        sym = np.searchsorted(vals, S) if alph else S
        est = all_estimators(sym, 2, include_predictors)
        return (float(min(est.values())) if est else float("nan")), est
    sym = np.searchsorted(vals, S)
    lit = all_estimators(sym, alph, include_predictors)
    if relabel_bits:
        word = max(1, (alph - 1).bit_length())
        bits = _to_bits(sym, word)
    else:
        # the bits of the samples as given, SP 800-90B section 3.1.3
        word = int(bits_per_symbol or max(1, int(S.max()).bit_length()))
        if int(S.max()) >> word:
            raise ValueError(f"a sample does not fit in {word} bits")
        bits = _to_bits(S, word)
    bit = all_estimators(bits, 2, include_predictors)
    h_original = float(min(lit.values()))
    h_bitstring = float(min(bit.values()))
    assessed = min(float(word), word * h_bitstring, h_original)
    per = {f"literal:{n}": v for n, v in lit.items()}
    per.update({f"bitstring:{n}": v for n, v in bit.items()})
    per.update(H_original=h_original, H_bitstring=h_bitstring, word_size=float(word))
    return float(assessed), per
