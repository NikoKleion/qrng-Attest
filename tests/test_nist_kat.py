# known-answer validation against the NIST SP 800-90B reference vectors and the output its ea_non_iid prints for them
import os
import re

import numpy as np

from _harness import skip
from qrng_attest import estimators as E
from qrng_attest import predictors as P

_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nist_vectors")
_NAMES = {"Most Common Value": "most_common_value", "Collision": "collision", "Markov": "markov",
          "Compression": "compression", "t-Tuple": "t_tuple", "LRS": "lrs", "MultiMCW Prediction": "multimcw",
          "Lag Prediction": "lag", "MultiMMC Prediction": "multimmc", "LZ78Y Prediction": "lz78y"}


def _load(name):
    path = os.path.join(_DIR, name)
    if not os.path.exists(path):
        skip(f"NIST vector {name} not present (download from the NIST repo to enable this KAT)")
    return np.frombuffer(open(path, "rb").read(), np.uint8).astype(int)


def _reference(name):
    # every min-entropy line of a reference output, keyed as min_entropy() keys its results
    path = os.path.join(_DIR, name)
    if not os.path.exists(path):
        skip(f"NIST reference output {name} not present")
    out, counts = {}, {}
    for line in open(path, encoding="utf-8"):
        m = re.match(r"\s*(Literal|Bitstring) (.+?) Estimate: min entropy = (\S+)", line)
        if m:
            out[f"{m.group(1).lower()}:{_NAMES[m.group(2)]}"] = float(m.group(3))
        m = re.match(r"\s*(Literal|Bitstring) (.+?) Estimate: (C|r|N) = (\d+)", line)
        if m:
            counts.setdefault(f"{m.group(1).lower()}:{_NAMES[m.group(2)]}", {})[m.group(3)] = int(m.group(4))
        m = re.match(r"\s*(H_original|H_bitstring|Assessed min entropy)\s*[:=]\s*(\S+)", line)
        if m:
            out["assessed" if m.group(1).startswith("Assessed") else m.group(1)] = float(m.group(2))
    return out, counts


def test_nist_mcv_micro_example():
    # MCV is hand-computable: L=100, 60 ones, H = -log2(p_u)
    S = np.array([1] * 60 + [0] * 40)
    assert abs(E.most_common_value(S, k=2) - 0.460167) < 2e-3


def test_nist_rand1_binary():
    # all ten estimators reproduce the reference tool to < 1e-9 bits
    ref, _ = _reference("rand1_short.res")
    est = E.all_estimators(_load("rand1_short.bin"), k=2)
    assert len(est) == 10
    for name, value in est.items():
        assert abs(value - ref[f"literal:{name}"]) < 1e-9, f"{name}: {value!r} vs {ref[f'literal:{name}']!r}"
    mn, _ = E.min_entropy(_load("rand1_short.bin"))
    assert abs(mn - ref["assessed"]) < 1e-9


def _multibit(vector, word):
    ref, _ = _reference(vector + ".res")
    mn, per = E.min_entropy(_load(vector + ".bin"))
    assert len(ref) == 6 + 4 + 3 + 4 + 3, sorted(ref)
    for key, value in ref.items():
        got = mn if key == "assessed" else per[key]
        assert abs(got - value) < 1e-9, f"{vector} {key}: {got!r} vs {value!r}"
    assert abs(mn - min(per["H_original"], word * per["H_bitstring"])) < 1e-9


def test_nist_rand4_multibit():
    # multi-bit assessment: literal (per-symbol) and bitstring (per-bit) tracks, every estimator of both
    _multibit("rand4_short", 4)


def test_nist_rand8_multibit():
    _multibit("rand8_short", 8)


def test_nist_predictor_counts():
    # correct predictions, longest run plus one and predictions made, as integers, for every predictor block
    for vector in ("rand1_short", "rand4_short", "rand8_short"):
        _, counts = _reference(vector + ".res")
        raw = _load(vector + ".bin")
        word = max(1, int(raw.max()).bit_length())
        bits = raw if word == 1 else E._to_bits(raw, word)
        checked = 0
        for key, want in counts.items():
            track, name = key.split(":")
            if name not in P.COUNTS:
                continue
            S, k = (raw, len(np.unique(raw))) if track == "literal" else (bits, 2)
            correct, n, run = P.COUNTS[name](S, k)
            assert (correct, run + 1, n) == (want["C"], want["r"], want["N"]), (vector, key, correct, run + 1, n, want)
            checked += 1
        assert checked == (4 if word == 1 else 8), (vector, checked)
