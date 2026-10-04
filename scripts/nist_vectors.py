# every estimator on the NIST reference vectors beside the reference output
#   python scripts/nist_vectors.py                     the three vendored vectors, writes results/nist_vectors.txt
#   python scripts/nist_vectors.py <clone> [vector]    vectors of a clone of usnistgov/SP800-90B_EntropyAssessment,
#                                                      one row each in results/nist_reference_vectors.txt
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from qrng_attest import estimators as E  # noqa: E402
from qrng_attest import predictors as P  # noqa: E402

VECTORS = os.path.join(ROOT, "tests", "nist_vectors")
NAMES = {"Most Common Value": "most_common_value", "Collision": "collision", "Markov": "markov",
         "Compression": "compression", "t-Tuple": "t_tuple", "LRS": "lrs", "MultiMCW Prediction": "multimcw",
         "Lag Prediction": "lag", "MultiMMC Prediction": "multimmc", "LZ78Y Prediction": "lz78y"}


def reference(name, folder=VECTORS):
    values, counts = {}, {}
    for line in open(os.path.join(folder, name), encoding="utf-8", errors="replace"):
        m = re.match(r"\s*(Literal|Bitstring) (.+?) Estimate: min entropy = (\S+)", line)
        if m:
            values[f"{m.group(1).lower()}:{NAMES[m.group(2)]}"] = float(m.group(3))
        m = re.match(r"\s*(Literal|Bitstring) (.+?) Estimate: (C|r|N) = (\d+)", line)
        if m:
            counts.setdefault(f"{m.group(1).lower()}:{NAMES[m.group(2)]}", {})[m.group(3)] = int(m.group(4))
        m = re.match(r"\s*(H_original|H_bitstring|Assessed min entropy)\s*[:=]\s*(\S+)", line)
        if m:
            values["assessed" if m.group(1).startswith("Assessed") else m.group(1)] = float(m.group(2))
    return values, counts


def main():
    out = ["Known-answer validation against the NIST SP 800-90B reference vectors",
           "reference values are the output of NIST ea_non_iid on the same vectors, tests/nist_vectors/*.res", ""]
    worst_all, counts_equal, blocks = 0.0, True, 0
    for vector, title in (("rand1_short", "binary"), ("rand4_short", "multi-bit, alphabet 16"),
                          ("rand8_short", "multi-bit, alphabet 256")):
        raw = np.frombuffer(open(os.path.join(VECTORS, vector + ".bin"), "rb").read(), np.uint8).astype(int)
        ref, counts = reference(vector + ".res")
        mn, per = E.min_entropy(raw)
        if title == "binary":
            per = {f"literal:{k}": v for k, v in per.items()}
            per["H_original"] = mn
        per["assessed"] = mn
        out.append(f"{vector} ({title})")
        out.append(f"  {'quantity':>30} {'this package':>14} {'NIST':>14} {'|difference|':>14}")
        worst = 0.0
        for key in sorted(ref):
            d = abs(per[key] - ref[key])
            worst = max(worst, d)
            out.append(f"  {key:>30} {per[key]:14.10f} {ref[key]:14.10f} {d:14.2e}")
        out.append(f"  worst deviation {worst:.2e}")
        word = max(1, int(raw.max()).bit_length())
        bits = raw if word == 1 else E._to_bits(raw, word)
        for key, want in sorted(counts.items()):
            track, name = key.split(":")
            if name in P.COUNTS:
                S, k = (raw, len(np.unique(raw))) if track == "literal" else (bits, 2)
                c, n, run = P.COUNTS[name](S, k)
                counts_equal &= (c, run + 1, n) == (want["C"], want["r"], want["N"])
                blocks += 1
        out.append("")
        worst_all = max(worst_all, worst)
    out.append(f"worst deviation over the three vectors {worst_all:.2e}")
    out.append(f"predictor counts C, r and N equal to the reference on {blocks} blocks: {counts_equal}")
    text = "\n".join(out) + "\n"
    with open(os.path.join(ROOT, "results", "nist_vectors.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    print(text, end="")


def clone_vectors(clone, names):
    # the full assessment on vectors of NIST's repository against cpp/selftest/refdata
    refdata = os.path.join(clone, "cpp", "selftest", "refdata")
    names = names or sorted(f[:-4] for f in os.listdir(refdata) if f.endswith(".res"))
    path = os.path.join(ROOT, "results", "nist_reference_vectors.txt")
    head = ["The SP 800-90B non-IID assessment on NIST's reference vectors, against the output NIST publishes for them",
            "every min entropy line of both tracks, H_original, H_bitstring and the assessed value",
            "",
            f"{'vector':22} {'samples':>9} {'alphabet':>9} {'values':>7} {'assessed':>14} {'NIST':>14} {'worst difference':>17}"]
    def saved():
        rows = {}
        if os.path.exists(path):
            for line in open(path, encoding="utf-8").read().splitlines()[len(head):]:
                if line.strip():
                    rows[line.split()[0]] = line
        return rows

    for name in names:
        raw = np.frombuffer(open(os.path.join(clone, "bin", name + ".bin"), "rb").read(), np.uint8).astype(int)
        ref, _ = reference(name + ".res", refdata)
        mn, per = E.min_entropy(raw)
        if "H_original" not in per:
            per = {f"literal:{k}": v for k, v in per.items()}
            per["H_original"] = mn
        per["assessed"] = mn
        worst = max(abs(per[key] - ref[key]) for key in ref)
        rows = saved()
        rows[name] = (f"{name:22} {len(raw):9d} {len(np.unique(raw)):9d} {len(ref):7d} {mn:14.10f} {ref['assessed']:14.10f} "
                      f"{worst:17.2e}")
        print(rows[name], flush=True)
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(head + [rows[k] for k in sorted(rows)]) + "\n")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        clone_vectors(sys.argv[1], sys.argv[2:])
    else:
        main()
