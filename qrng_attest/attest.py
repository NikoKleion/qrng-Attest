# device attestation: classical and readout-corrected estimates, correlation and drift
import math
import numpy as np

from .estimators import min_entropy, _upper_bound_p
from .quantum import (QuantumSource, QuantumArraySource, attested_qubit_entropy, attested_coherence_entropy,
                      shot_entropy, stokes_lower, stokes_upper, worst_entropy)
from .health import health_report


# cross-qubit geometry
def mutual_information(x, y):
    # mutual information in bits between two binary streams.
    x = np.asarray(x, int); y = np.asarray(y, int); n = len(x)
    I = 0.0
    for a in (0, 1):
        px = np.mean(x == a)
        for b in (0, 1):
            py = np.mean(y == b); pxy = np.mean((x == a) & (y == b))
            if pxy > 0 and px > 0 and py > 0:
                I += pxy * math.log2(pxy / (px * py))
    return max(I, 0.0)


def independence(shots, threshold=0.005):
    # largest and mean mutual information over all qubit pairs. Change threshold for the independence cutoff.
    shots = np.asarray(shots, int); nq = shots.shape[1]
    pairs = []
    for i in range(nq):
        for j in range(i + 1, nq):
            pairs.append(((i, j), mutual_information(shots[:, i], shots[:, j])))
    max_mi = max((mi for _, mi in pairs), default=0.0)
    return {"max_mutual_info": max_mi, "mean_mutual_info": float(np.mean([mi for _, mi in pairs])) if pairs else 0.0,
            "independent": max_mi < threshold,
            "worst_pair": max(pairs, key=lambda kv: kv[1])[0] if pairs else None}


def pairwise_joint_bound(shots):
    # crosstalk penalty: a per-bit bound on the joint entropy rate from the worst qubit pair.
    shots = np.asarray(shots, int)
    n, nq = shots.shape
    if nq < 2 or n < 2:
        return 1.0, None
    worst = 1.0; worst_pair = None
    for i in range(nq):
        for j in range(i + 1, nq):
            joint = shots[:, i] * 2 + shots[:, j]
            phat = float(np.bincount(joint, minlength=4).max()) / n
            h_pair = -math.log2(_upper_bound_p(phat, n))
            rate = (h_pair + (nq - 2)) / nq
            if rate < worst:
                worst = rate; worst_pair = (i, j)
    return float(worst), worst_pair


def drift(bits, windows=20):
    # drift check: compare the spread of per-window bias to the spread expected from sampling alone.
    bits = np.asarray(bits, int); n = len(bits)
    w = n // windows
    if w < 30:
        return {"drift_score": 1.0, "trend_slope": 0.0, "stationary": True}
    p0 = np.array([np.mean(bits[k * w:(k + 1) * w] == 0) for k in range(windows)])
    pbar = float(np.mean(p0))
    expected_var = pbar * (1 - pbar) / w
    observed_var = float(np.var(p0))
    score = observed_var / expected_var if expected_var > 0 else 1.0
    slope = float(np.polyfit(np.arange(windows), p0, 1)[0])
    return {"drift_score": score, "trend_slope": slope, "stationary": score < 3.0}


# the attestation
class QuantumAttestation:
    def __init__(self, per_qubit, indep, drift_res, describe, attested):
        self.per_qubit = per_qubit
        self.independence = indep
        self.drift = drift_res
        self.describe = describe
        self.attested_min_entropy = attested
        # without X and Y shots the value holds for a pure source state only, so the source is not passed
        self.attested = all(q["coherence_low"] is not None for q in per_qubit)
        self.passed = (attested > 0.01 and self.attested
                       and all(q["health"]["passed"] for q in per_qubit)
                       and indep["independent"] and drift_res["stationary"]
                       and not any(q["readout_dominated"] for q in per_qubit))

    def shot_bound(self, epsilon):
        # smooth min-entropy of the whole output: (bits, alpha, the term that sets it), capped by the classical values
        if any(q["coherence_low"] is None for q in self.per_qubit):
            return None
        n = self.describe["n_shots"]
        bits, order = shot_entropy([(q["coherence_low"], q["s3_low"], q["s3_high"]) for q in self.per_qubit],
                                   n, epsilon)
        # each qubit's SP 800-90B value caps its own stream; the joint bound caps the correlation between them
        terms = {"the alpha entropy": bits, "the SP 800-90B values": n * sum(q["h_classical"] for q in self.per_qubit),
                 "the bound for correlated qubits": n * len(self.per_qubit) * self.independence["joint_rate_bound"]}
        term = min(terms, key=terms.get)
        return terms[term], order, term

    def extraction_plan(self, epsilon=2.0 ** -64):
        from .extract import plan
        n = self.describe["n_shots"] * self.describe["n_qubits"]
        one = plan(n, self.attested_min_entropy, epsilon=epsilon,
                   note="rate is the attested per-bit bound, not the observed one")
        many = self.shot_bound(epsilon / 2.0)
        if many is None or many[1] is None:
            return one
        ell = math.floor(many[0] - 2.0 * math.log2(1.0 / epsilon))
        if ell <= one.output_bits:
            return one
        out = plan(n, self.attested_min_entropy, epsilon=epsilon, ell=ell, order=many[1],
                   note="output length from the shots taken together")
        out.term = many[2]
        return out

    def __str__(self):
        d = self.describe
        L = [f"Quantum entropy attestation: {d['name']}  [{d['provenance']}]",
             f"  {d['n_shots']} shots x {d['n_qubits']} qubits, expected source state |{d['expected_state']}>",
             f"  attested min-entropy = {self.attested_min_entropy:.4f} bits/bit   "
             f"(verdict: {'PASS' if self.passed else 'FLAGGED'})",
             "  per qubit (classical non-IID vs readout-corrected quantum):"]
        for i, q in enumerate(self.per_qubit):
            tag = "  READOUT-DOMINATED" if q["readout_dominated"] else ""
            hh = "ok" if q["health"]["passed"] else "HEALTH-FAIL"
            coh = "" if q["coherence"] is None else f"  coherence {q['coherence']:.3f} (at least {q['coherence_low']:.3f})"
            L.append(f"    q{i}: classical {q['h_classical']:.3f}  quantum {q['h_quantum']:.3f}  "
                     f"r_z={q['r_z']:+.3f}{coh}  [{hh}]{tag}")
        im = self.independence
        L.append(f"  independence: max pairwise mutual info {im['max_mutual_info']:.4f} bits "
                 f"({'independent' if im['independent'] else 'CORRELATED, crosstalk'})")
        if im.get("joint_worst_pair") is not None:
            L.append(f"  joint rate bound: {im['joint_rate_bound']:.4f} bits/bit "
                     f"(worst pair q{im['joint_worst_pair'][0]}, q{im['joint_worst_pair'][1]})")
        dr = self.drift
        L.append(f"  stationarity: drift score {dr['drift_score']:.1f} "
                 f"({'stationary' if dr['stationary'] else 'DRIFTING'}), trend slope {dr['trend_slope']:+.2e}/window")
        for f in d["flags"]:
            L.append(f"  note: {f}")
        pl = self.extraction_plan()
        L.append(f"  extractable: {pl.output_bits} uniform bits from {pl.n_input} attested bits "
                 f"(epsilon 2^-64); a 256-bit key needs {pl.bits_for(256)} bits")
        for r in remediation(self):
            L.append(f"  action: {r}")
        for r in advice(self):
            L.append(f"  advice: {r}")
        L.append("  assumes a trusted measurement and a characterized readout"
                 + ("" if all(q["coherence"] is not None for q in self.per_qubit) else ", and a pure source state"))
        return "\n".join(L)


def qubit_subset(source, include_predictors=False, min_qubits=1, run_iid=False):
    # drop the worst qubit one at a time and report which subset gives the most total bits per shot.
    src = source if isinstance(source, QuantumSource) else QuantumArraySource(np.asarray(source))
    shots = np.asarray(src.shots(), int)
    Ms = src.readout_matrices()
    cal = getattr(src, "calibration_shots", None)
    basis = src.basis_shots()
    part = lambda key, idx: None if basis is None or basis.get(key) is None else basis[key][:, idx]
    # a qubit's classical entropy does not depend on the others, so compute it once and reuse it
    ch_full = {q: _classical_entropy(shots[:, q], include_predictors, run_iid) for q in range(shots.shape[1])}

    def _att(idx):
        return attest(QuantumArraySource(shots[:, idx], readout=None if Ms is None else [Ms[q] for q in idx],
                                         calibration_shots=cal, x_shots=part("x", idx), y_shots=part("y", idx)),
                      include_predictors=include_predictors, classical_h=[ch_full[q] for q in idx])

    keep = list(range(shots.shape[1]))
    base = _att(keep)
    steps = [{"qubits": list(keep), "rate": base.attested_min_entropy, "dropped": None}]
    while len(keep) > min_qubits:
        best = None
        for q in keep:
            trial = [x for x in keep if x != q]
            a = _att(trial)
            if best is None or a.attested_min_entropy > best[1]:
                best = (q, a.attested_min_entropy, trial)
        if best is None or best[1] <= steps[-1]["rate"] + 1e-12:
            break
        keep = best[2]
        steps.append({"qubits": list(keep), "rate": best[1], "dropped": best[0]})
    for st in steps:
        st["bits_per_shot"] = st["rate"] * len(st["qubits"])
    best = max(steps, key=lambda st: st["bits_per_shot"])
    base_total = steps[0]["bits_per_shot"]
    return {"best_qubits": best["qubits"], "best_rate": best["rate"],
            "base_rate": steps[0]["rate"], "steps": steps,
            "dropped": [q for q in steps[0]["qubits"] if q not in best["qubits"]],
            "total_bits_per_shot": best["bits_per_shot"], "base_bits_per_shot": base_total,
            "worth_it": best["bits_per_shot"] > base_total + 1e-12}


def remediation(att):
    # the failed checks as a list of actions, or one line saying there is none
    problems = []
    v = att.attested_min_entropy
    tol = 1e-6
    if any(q["coherence_low"] is None for q in att.per_qubit):
        problems.append("no X or Y basis shots: the value holds for a pure source state only; take them for an "
                        "attested bound")
    for i, q in enumerate(att.per_qubit):
        if q["readout_dominated"]:
            if q["h_population"] <= v + tol and q["h_population"] <= q["h_classical"] + tol:
                problems.append(f"q{i}: readout error sets the attested rate ({q['h_population']:.2f} bits/bit); "
                                f"recalibrate the discriminator or exclude this qubit")
            else:
                gap = q["h_observed"] - q["h_population"]
                problems.append(f"q{i}: readout error costs {gap:.2f} bits on this qubit but is not the current "
                                f"limit; recalibrate it before you raise the shot count")
        if q["h_coherence"] is not None and q["h_coherence"] <= 0.01:
            problems.append(f"q{i}: no coherence is certified (at least {q['coherence_low']:.3f}); check the state "
                            f"preparation and the X and Y basis circuits")
        if not q["health"]["passed"]:
            problems.append(f"q{i}: health test failed; discard recent output and run again after recalibration")
    im = att.independence
    if not im["independent"]:
        p = im.get("joint_worst_pair") or im.get("worst_pair")
        if p is not None:
            problems.append(f"q{p[0]} and q{p[1]} are correlated; separate them physically or drop one from the pool")
    if not att.drift["stationary"]:
        problems.append("source is drifting; shorten the re-attestation interval and re-run calibration")

    return problems or ["no action required; re-attest after the next calibration cycle"]


def advice(att):
    # what sets the attested rate and how to raise it; none of these is a failed check
    out = []
    v = att.attested_min_entropy
    tol = 1e-6
    for i, q in enumerate(att.per_qubit):
        if (q["h_coherence"] is not None and 0.01 < q["h_coherence"] <= v + tol
                and q["h_coherence"] < q["h_population"] - tol):
            out.append(f"q{i}: the coherence bound sets the attested rate ({q['h_coherence']:.2f} bits/bit at "
                       f"coherence {q['coherence_low']:.3f}); raise the X and Y basis and calibration shot "
                       f"counts, or check the state preparation")
    rates = [q["h_quantum"] for q in att.per_qubit]
    if len(rates) > 1:
        worst = int(np.argmin(rates))
        if rates[worst] < 0.5 * float(np.median(rates)):
            out.append(f"q{worst} yields {rates[worst]:.2f} bits/bit against a median of "
                       f"{float(np.median(rates)):.2f}; check whether dropping it raises total throughput")
    # the statistical estimator, not the device, set the bound
    min_hq = min((q["h_quantum"] for q in att.per_qubit), default=1.0)
    joint = att.independence.get("joint_rate_bound", 1.0)
    if v < min_hq - tol and v < joint - tol:
        n = att.describe.get("n_shots", 0)
        out.append(f"attested rate {v:.2f} is set by the statistical estimator at {n} shots, not the device; "
                   f"increase shots to raise it")
    return out


# IID track budget for the per-qubit classical estimate.
IID_PERMS = 10000
IID_MAX = 8000


def _classical_entropy(stream, include_predictors, run_iid):
    # SP 800-90B classical min-entropy for one qubit stream, returned as (bits, track).
    h_noniid, _ = min_entropy(stream, k=2, include_predictors=include_predictors)
    if not run_iid:
        return h_noniid, "non-IID"
    from .iid import chi_square_tests, lrs_test, permutation_test, _mcv_assessed
    # skip the permutation test when the fast checks already reject IID
    chi_ok, _pi, _pg = chi_square_tests(stream)
    lrs_ok, _pl = lrs_test(stream)
    if not (chi_ok and lrs_ok):
        return h_noniid, "non-IID"
    if not permutation_test(stream, perms=IID_PERMS, max_samples=IID_MAX).is_iid:
        return h_noniid, "non-IID"
    assessed, _ = _mcv_assessed(stream)
    return assessed, "IID"


def attest(source, include_predictors=True, run_iid=False, classical_h=None):
    # classical_h: per-qubit (bits, track), computed once for the qubit-subset search
    src = source if isinstance(source, QuantumSource) else QuantumArraySource(np.asarray(source))
    shots = np.asarray(src.shots(), int)
    nshots, nq = shots.shape
    Ms = src.readout_matrices()
    basis = src.basis_shots()
    cal = getattr(src, "calibration_shots", None)
    per_qubit = []
    attested = math.inf
    for q in range(nq):
        stream = shots[:, q]
        if classical_h is not None:
            h_classical, track = classical_h[q]
        else:
            h_classical, track = _classical_entropy(stream, include_predictors, run_iid)
        qe = attested_qubit_entropy(stream, Ms[q] if Ms is not None else None, calib_shots=cal)
        ce = None
        if basis is not None:
            ce = attested_coherence_entropy(basis["x"][:, q], None if basis["y"] is None else basis["y"][:, q],
                                            Ms[q] if Ms is not None else None, calib_shots=cal)
        h_quantum = qe["h_quantum"] if ce is None else min(qe["h_quantum"], ce["h_coherence"])
        Mq = Ms[q] if Ms is not None else None
        s3 = (stokes_lower(stream, Mq, calib_shots=cal)[1], stokes_upper(stream, Mq, calib_shots=cal))
        h_q = min(h_classical, h_quantum)
        hr = health_report(stream, max(h_q, 1e-3))
        per_qubit.append({"h_classical": h_classical, "h_quantum": h_quantum, "r_z": qe["r_z"],
                          "h_observed": qe["h_observed"], "classical_track": track,
                          "readout_dominated": qe["readout_dominated"], "health": hr,
                          "h_population": qe["h_quantum"],
                          "h_coherence": None if ce is None else ce["h_coherence"],
                          "s3_low": s3[0], "s3_high": s3[1],
                          "h_long_run": None if ce is None else worst_entropy(ce["coherence_low"], s3[0], s3[1], 1.0),
                          "coherence": None if ce is None else ce["coherence"],
                          "coherence_low": None if ce is None else ce["coherence_low"]})
        attested = min(attested, h_q)
    indep = independence(shots)
    drift_res = drift(shots[:, 0])
    joint_rate, joint_pair = pairwise_joint_bound(shots)
    indep["joint_rate_bound"] = joint_rate
    indep["joint_worst_pair"] = joint_pair
    attested = min(attested, joint_rate)
    return QuantumAttestation(per_qubit, indep, drift_res, src.describe(), float(attested))
