# the four SP 800-90B section 6.3 predictors, a translation of NIST's reference implementation
# (usnistgov/SP800-90B_EntropyAssessment 1.1.8: cpp/non_iid/*.h and cpp/shared/utils.h); see NOTICE
import math

import numpy as np

ZALPHA = 2.5758293035489008          # utils.h ZALPHA
ITERMAX = 1076                       # utils.h ITERMAX
D_LAG = 128                          # lag_test.h D_LAG
MCW_WINDOWS = (63, 255, 1023, 4095)  # multi_mcw_test.h W
D_MMC = 16                           # multi_mmc_test.h D_MMC
MMC_MAX_ENTRIES = 100000             # multi_mmc_test.h MAX_ENTRIES
B_LZ = 16                            # lz78y_test.h B_len
LZ_MAX_DICTIONARY = 65536            # lz78y_test.h MAX_DICTIONARY_SIZE
BLOCK = 1 << 19                      # cells of a score table handled at once by the array forms


def _run_function(p, r, n):
    # utils.h prediction_estimate_function: log of the chance of no run of r correct predictions in n trials
    q = 1.0 - p
    x, last = 1.0, 0.0
    i = 0
    while i <= 65 and (x - last) > np.finfo(float).eps * x:
        last = x
        x = 1.0 + q * math.exp(r * math.log(p) + (r + 1.0) * math.log(x))
        i += 1
    a, b = 1.0 - p * x, (r + 1.0 - r * x) * q
    if a <= 0.0 or b <= 0.0:
        return -math.inf if a <= 0.0 else math.inf
    return math.log(a) - math.log(b) - (n + 1.0) * math.log(x)


def _close(a, b):
    # utils.h relEpsilonEqual at 4 units in the last place
    if a == b:
        return True
    return abs(a - b) <= 4 * np.spacing(max(abs(a), abs(b)))


def p_local(max_run, n, low):
    # utils.h calc_p_local: the rate at which the longest run would be max_run + 1 with probability 0.01
    target = math.log(0.99)
    lo, hi = low, 1.0
    lo_val, hi_val = math.inf, -math.inf
    p = (lo + hi) / 2.0
    val = _run_function(p, max_run + 1, n)
    for _ in range(ITERMAX):
        if _close(val, target):
            break
        if target < val:
            lo, lo_val = p, val
        else:
            hi, hi_val = p, val
        if lo >= hi:
            p = min(max(lo, hi), 1.0)
            break
        if not (low <= lo <= 1.0 and low <= hi <= 1.0):
            p = 1.0
            break
        if not (min(lo_val, hi_val) <= target <= max(lo_val, hi_val)):
            p = 1.0
            break
        last = p
        p = (lo + hi) / 2.0
        if not (lo < p < hi) or last == p:
            p = hi
            break
        val = _run_function(p, max_run + 1, n)
        if not (min(lo_val, hi_val) <= val <= max(lo_val, hi_val)):
            p = hi
            break
    return p


def prediction_estimate(correct, n, max_run, k):
    # utils.h predictionEstimate: min-entropy from the count of correct predictions and the longest run of them
    cur = 1.0 / k
    pg = correct / n
    if pg > 0:
        bound = min(1.0, pg + ZALPHA * math.sqrt(pg * (1.0 - pg) / (n - 1.0)))
    else:
        bound = 1.0 - 0.01 ** (1.0 / n)
    cur = max(cur, bound)
    if cur < 1.0 and _run_function(cur, max_run + 1, n) > math.log(0.99):
        cur = max(cur, p_local(max_run, n, cur))
    return -math.log2(cur)


def lag_loop(S):
    # lag_test.h, one sample at a time: (correct predictions, predictions made, longest run of correct ones)
    S = np.asarray(S)
    n = len(S)
    if n < 3:
        return None
    scores = np.zeros(D_LAG, np.int64)
    winner, high = 0, 0
    correct = run = best = 0
    for i in range(1, n):
        cur = S[i]
        if cur == S[i - winner - 1]:
            correct += 1
            run += 1
            if run > best:
                best = run
        else:
            run = 0
        hit = np.flatnonzero(S[max(0, i - D_LAG):i][::-1] == cur)   # lags, smallest first
        if hit.size:
            scores[hit] += 1
            top = scores[hit]
            m = top.max()
            if m >= high:
                high = int(m)
                winner = int(hit[top == m][-1])
    return correct, n - 1, best


def multimcw_loop(S, k):
    # multi_mcw_test.h, one sample at a time
    data = [int(x) for x in S]
    n = len(data)
    W = MCW_WINDOWS
    if n < W[-1] + 1:
        return None
    cnt = [[0] * k for _ in W]
    pos = [[0] * k for _ in W]
    top = [0] * len(W)
    freq = [0] * len(W)
    score = [0] * len(W)
    for i in range(W[-1]):
        x = data[i]
        for j in range(len(W)):
            if i < W[j]:
                cnt[j][x] += 1
                if top[j] <= cnt[j][x]:
                    top[j] = cnt[j][x]
                    freq[j] = x
                pos[j][x] = i
    winner = 0
    correct = run = best = 0
    for i in range(W[0], n):
        x = data[i]
        if freq[winner] == x:
            correct += 1
            run += 1
            if run > best:
                best = run
        else:
            run = 0
        for j in range(len(W)):
            if i >= W[j] and freq[j] == x:
                score[j] += 1
                if score[j] >= score[winner]:
                    winner = j
        for j in range(len(W)):
            if i >= W[j]:
                out = data[i - W[j]]
                cj, pj = cnt[j], pos[j]
                cj[out] -= 1
                cj[x] += 1
                pj[x] = i
                if out != freq[j]:
                    if top[j] <= cj[x]:
                        top[j] = cj[x]
                        freq[j] = x
                else:
                    m, f, mp = top[j] - 1, freq[j], i - W[j]
                    for s in range(k):
                        c = cj[s]
                        if m < c or (m == c and mp <= pj[s]):
                            m, f, mp = c, s, pj[s]
                    top[j], freq[j] = m, f
    return correct, n - W[0], best


def _bump(entry, y, make_new):
    # utils.h PostfixDictionary::incrementPostfix; entry is [counts, best count, prediction]
    counts = entry[0]
    c = counts.get(y)
    if c is not None:
        c += 1
        counts[y] = c
        new = False
    elif make_new:
        c = 1
        counts[y] = 1
        new = True
    else:
        return False
    if c > entry[1] or (c == entry[1] and y > entry[2]):
        entry[1], entry[2] = c, y
    return new


def multimmc_loop(S, k):
    # multi_mmc_test.h, one sample at a time
    data = [int(x) for x in S]
    n = len(data)
    if n < 3:
        return None
    N = n - 2
    power = [k ** d for d in range(D_MMC)]
    M = [dict() for _ in range(D_MMC)]
    entries = [0] * D_MMC
    score = [0] * D_MMC
    code = 0
    for d in range(D_MMC):
        if d < N:
            code = code * k + data[d]          # data[0..d], the latest symbol in the lowest place
            _bump(M[d].setdefault(code, [{}, 0, 0]), data[d + 1], True)
            entries[d] = 1
    winner = 0
    correct = run = best = 0
    for i in range(2, n):
        y = data[i]
        found = False
        cur_winner = winner
        code = 0
        entry = None
        for d in range(min(D_MMC, i - 1)):
            code += data[i - d - 1] * power[d]   # data[i-d-1..i-1]
            if d == 0 or found:
                entry = M[d].get(code)
                found = entry is not None
            if found:
                if entry[2] == y:
                    score[d] += 1
                    if score[d] >= score[winner]:
                        winner = d
                    if d == cur_winner:
                        correct += 1
                        run += 1
                        if run > best:
                            best = run
                elif d == cur_winner:
                    run = 0
                if _bump(entry, y, entries[d] < MMC_MAX_ENTRIES):
                    entries[d] += 1
            elif entries[d] < MMC_MAX_ENTRIES:
                _bump(M[d].setdefault(code, [{}, 0, 0]), y, True)
                entries[d] += 1
    return correct, N, best


def lz78y_loop(S, k):
    # lz78y_test.h, one sample at a time
    data = [int(x) for x in S]
    n = len(data)
    if n < B_LZ + 2:
        return None
    N = n - B_LZ - 1
    power = [k ** j for j in range(B_LZ)]
    D = [dict() for _ in range(B_LZ)]
    size = 0
    code = 0
    for j in range(1, B_LZ + 1):
        code += data[B_LZ - j] * power[j - 1]      # data[B-j..B-1]
        _bump(D[j - 1].setdefault(code, [{}, 0, 0]), data[B_LZ], True)
        size += 1
    correct = run = best = 0
    codes = [0] * (B_LZ + 1)
    for i in range(B_LZ + 1, n):
        y = data[i]
        code = 0
        for j in range(1, B_LZ + 1):
            code += data[i - j] * power[j - 1]     # data[i-j..i-1]
            codes[j] = code
        have, guess, most = False, 0, 0
        for j in range(B_LZ, 0, -1):
            entry = D[j - 1].get(codes[j])
            if entry is not None:
                if entry[1] > most:
                    most, guess, have = entry[1], entry[2], True
                _bump(entry, y, True)
            elif size < LZ_MAX_DICTIONARY:
                _bump(D[j - 1].setdefault(codes[j], [{}, 0, 0]), y, True)
                size += 1
        if have and guess == y:
            correct += 1
            run += 1
            if run > best:
                best = run
        else:
            run = 0
    return correct, N, best


def _longest_run(hit):
    # length of the longest run of True
    edge = np.flatnonzero(np.diff(np.concatenate(([0], np.asarray(hit, np.int8), [0]))))
    return int((edge[1::2] - edge[::2]).max()) if edge.size else 0


def _leaders(hits, scores, leader):
    # leader before each step, then scores and leader after; a step's new leader is its last hit with the top score
    rows, steps = hits.shape
    after = scores[:, None] + np.cumsum(hits, axis=1, dtype=np.int32)
    best = hits & (after == after.max(axis=0))
    last = rows - 1 - np.argmax(best[::-1], axis=0)
    at = np.maximum.accumulate(np.where(best[last, np.arange(steps)], np.arange(steps), -1))
    led = np.where(at >= 0, last[np.maximum(at, 0)], leader)
    return np.concatenate(([leader], led[:-1])), after[:, -1].copy(), int(led[-1])


def _earlier(keys, one, zeros, ones, keep=None):
    # how often 0 and 1 followed each position's key before it; the two tables carry the counts between blocks
    m = len(keys)
    order = np.argsort(keys, kind="stable")
    ks = keys[order]
    idx = np.arange(m)
    start = np.maximum.accumulate(np.where(np.concatenate(([True], ks[1:] != ks[:-1])), idx, 0))
    w1 = one if keep is None else one & keep
    w0 = ~one if keep is None else ~one & keep
    run = np.cumsum(w1[order]) - w1[order]
    n1 = np.empty(m, np.int64)
    n1[order] = run - run[start]
    n0 = np.empty(m, np.int64)
    if keep is None:
        n0[order] = idx - start
        n0 -= n1
    else:
        run = np.cumsum(w0[order]) - w0[order]
        n0[order] = run - run[start]
    n0 += zeros[keys]
    n1 += ones[keys]
    zeros += np.bincount(keys[w0], minlength=len(zeros))
    ones += np.bincount(keys[w1], minlength=len(ones))
    return n0, n1


def lag_counts(S):
    # lag_test.h in blocks of steps: (correct predictions, predictions made, longest run of correct ones)
    S = np.asarray(S)
    n = len(S)
    if n < 3:
        return None
    S = S.astype(np.int8 if 0 <= S.min() and S.max() < 127 else np.int64)
    past = np.lib.stride_tricks.sliding_window_view(np.concatenate((np.full(D_LAG, S.min() - 1, S.dtype), S)), D_LAG)
    scores, leader = np.zeros(D_LAG, np.int32), 0
    correct = np.zeros(n, bool)
    step = max(1, BLOCK // D_LAG)
    for a in range(1, n, step):
        b = min(n, a + step)
        hits = np.ascontiguousarray((past[a:b, ::-1] == S[a:b, None]).T)     # row d: S[t] == S[t - d - 1]
        prior, scores, leader = _leaders(hits, scores, leader)
        correct[a:b] = hits[prior, np.arange(b - a)]
    return int(correct.sum()), n - 1, _longest_run(correct)


def _multimcw_bits(S):
    # multi_mcw_test.h on bits: the most common value of a window of odd length is its majority
    n = len(S)
    W = MCW_WINDOWS
    if n < W[-1] + 1:
        return None
    one = S.astype(bool)
    scores, leader = np.zeros(len(W), np.int32), 0
    correct = np.zeros(n, bool)
    step = max(1, BLOCK // len(W))
    for a in range(W[0], n, step):
        b = min(n, a + step)
        lo = max(0, a - W[-1])
        ones = np.concatenate(([0], np.cumsum(one[lo:b])))                   # ones[t - lo]: ones in S[lo:t]
        hits = np.zeros((len(W), b - a), bool)
        for j, w in enumerate(W):
            t = np.arange(max(a, w), b)
            hits[j, t - a] = (2 * (ones[t - lo] - ones[t - w - lo]) > w) == one[t]
        prior, scores, leader = _leaders(hits, scores, leader)
        correct[a:b] = hits[prior, np.arange(b - a)]
    return int(correct.sum()), n - W[0], _longest_run(correct)


def _multimmc_bits(S):
    # multi_mmc_test.h on bits: the counts after each context by sorting; the entry limit binds at the last depth
    n = len(S)
    if n < 3:
        return None
    depths = min(D_MMC, n - 2)
    bit, one = S.astype(np.uint16), S.astype(bool)
    counts = [(np.zeros(2 << d, np.int64), np.zeros(2 << d, np.int64)) for d in range(depths)]
    limited = D_MMC - 1 if (2 << D_MMC) > MMC_MAX_ENTRIES else -1
    entered, room = np.zeros(2 << D_MMC, bool), MMC_MAX_ENTRIES
    scores, leader = np.zeros(depths, np.int32), 0
    status = np.zeros(n, np.uint8)                       # per step: 0 the leader's context is not held, 1 wrong, 2 correct
    step = max(1, BLOCK // depths)
    for a in range(1, n, step):
        b = min(n, a + step)
        state = np.zeros((depths, b - a), np.uint8)
        ctx = np.zeros(b - a, np.uint16)
        for d in range(depths):
            lo = max(a, d + 1)                           # depth d reads S[t - d - 1..t - 1]
            if lo >= b:
                break
            ctx[lo - a:] |= bit[lo - d - 1:b - d - 1] << d
            keys, y, keep = ctx[lo - a:], one[lo:b], None
            if d == limited:
                pairs = 2 * keys.astype(np.int64) + y
                if room:                                 # pairs enter in order of first appearance
                    new, first = np.unique(pairs, return_index=True)
                    fresh = ~entered[new]
                    take = new[fresh][np.argsort(first[fresh])[:room]]
                    entered[take] = True
                    room -= len(take)
                keep = entered[pairs]
            n0, n1 = _earlier(keys, y, *counts[d], keep)
            state[d, lo - a:] = np.where(n0 + n1 > 0, 1 + ((n1 >= n0) == y), 0)
        prior, scores, leader = _leaders(state == 2, scores, leader)
        status[a:b] = state[prior, np.arange(b - a)]
    status = status[2:]
    correct = status == 2
    return int(correct.sum()), n - 2, _longest_run(correct[status > 0])


def _lz78y_bits(S):
    # lz78y_test.h on bits: contexts enter the dictionary in order of time, longest first, up to its limit
    n = len(S)
    if n < B_LZ + 2:
        return None
    bit, one = S.astype(np.uint16), S.astype(bool)
    lengths = range(1, B_LZ + 1)
    counts = [(np.zeros(1 << j, np.int64), np.zeros(1 << j, np.int64)) for j in lengths]
    held = [np.zeros(1 << j, bool) for j in lengths]
    code = 0
    for j in lengths:                                    # the entries made before the first prediction
        code |= int(bit[B_LZ - j]) << (j - 1)
        held[j - 1][code] = True
    room = max(0, LZ_MAX_DICTIONARY - B_LZ)
    correct = np.zeros(n, bool)
    step = max(1, BLOCK // B_LZ)
    for a in range(B_LZ, n, step):
        b = min(n, a + step)
        y = one[a:b]
        if room:
            ctx, new = np.zeros(b - a, np.uint16), []
            for j in lengths:
                ctx |= bit[a - j:b - j] << (j - 1)
                codes, first = np.unique(ctx, return_index=True)
                fresh = ~held[j - 1][codes]
                new.append((codes[fresh], first[fresh] * (B_LZ + 1) + B_LZ - j))
            rank = np.concatenate([r for _, r in new])
            limit = np.partition(rank, room - 1)[room - 1] if len(rank) > room else np.inf
            for j, (codes, r) in zip(lengths, new):
                held[j - 1][codes[r <= limit]] = True
            room -= min(room, len(rank))
        most, guess = np.zeros(b - a, np.int64), np.zeros(b - a, bool)
        ctx = np.zeros(b - a, np.uint16)
        for j in lengths:                                # the largest count wins, the longest context at equal counts
            ctx |= bit[a - j:b - j] << (j - 1)
            n0, n1 = _earlier(ctx, y, *counts[j - 1])
            best = np.maximum(n0, n1)
            use = held[j - 1][ctx] & (best > 0) & (best >= most)
            most[use], guess[use] = best[use], (n1 >= n0)[use]
        correct[a:b] = (most > 0) & (guess == y)
    correct = correct[B_LZ + 1:]
    return int(correct.sum()), n - B_LZ - 1, _longest_run(correct)


def _are_bits(S, k):
    S = np.asarray(S)
    return k == 2 and S.size > 0 and S.min() >= 0 and S.max() <= 1


def multimcw_counts(S, k):
    # the array form on bits where every window is odd, the loop otherwise
    if _are_bits(S, k) and all(w % 2 for w in MCW_WINDOWS):
        return _multimcw_bits(np.asarray(S))
    return multimcw_loop(S, k)


def multimmc_counts(S, k):
    # the array form on bits where the entry limit cannot bind above the last depth, the loop otherwise
    if _are_bits(S, k) and D_MMC <= 16 and 2 ** D_MMC <= MMC_MAX_ENTRIES:
        return _multimmc_bits(np.asarray(S))
    return multimmc_loop(S, k)


def lz78y_counts(S, k):
    # the array form on bits, the loop otherwise
    if _are_bits(S, k) and B_LZ <= 16:
        return _lz78y_bits(np.asarray(S))
    return lz78y_loop(S, k)


COUNTS = {"multimcw": multimcw_counts, "lag": lambda S, k: lag_counts(S), "multimmc": multimmc_counts,
          "lz78y": lz78y_counts}
LOOPS = {"multimcw": multimcw_loop, "lag": lambda S, k: lag_loop(S), "multimmc": multimmc_loop, "lz78y": lz78y_loop}


def estimate(name, S, k):
    # min-entropy of one predictor, or None where the reference skips it for lack of samples
    counts = COUNTS[name](S, k)
    if counts is None or counts[1] < 2:
        return None
    return prediction_estimate(counts[0], counts[1], counts[2], k)
