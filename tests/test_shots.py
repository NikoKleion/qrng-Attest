# the bound on n shots taken together: the alpha entropy of a measured qubit, its worst case, and the output length
import math

import numpy as np

import qrng_attest as es
from qrng_attest import quantum as Q
from qrng_attest.attest import attest

M = Q.assignment_matrix(0.99, 0.98)


def _powm(A, x):
    w, v = np.linalg.eigh(A)
    keep = w > 1e-13 * w.max()
    return (v[:, keep] * w[keep] ** x) @ v[:, keep].conj().T


def _by_definition(s1, s2, s3, alpha):
    # Tomamichel, Colbeck and Renner 2009, eq. 16, by matrix powers on a purification of the qubit
    rho = 0.5 * np.array([[1 + s3, s1 - 1j * s2], [s1 + 1j * s2, 1 - s3]])
    psi = np.kron(_powm(rho, 0.5), np.eye(2)) @ np.array([1, 0, 0, 1], complex)
    state = np.zeros((4, 4), complex)
    for z in (0, 1):
        v = np.kron(np.diag([1.0 - z, float(z)]), np.eye(2)) @ psi
        state += np.outer(v, v.conj())
    ra = _powm(state, alpha)
    sigma = _powm(ra.reshape(2, 2, 2, 2).trace(axis1=0, axis2=2), 1.0 / alpha)
    sigma = sigma / np.trace(sigma).real
    return float(np.log2(np.trace(ra @ np.kron(np.eye(2), _powm(sigma, 1.0 - alpha))).real) / (1.0 - alpha))


def test_alpha_entropy_equals_its_definition():
    rng = np.random.default_rng(1)
    for _ in range(60):
        v = rng.normal(size=3)
        s1, s2, s3 = 0.999 * rng.random() ** (1 / 3) * v / np.linalg.norm(v)
        alpha = 1.001 + 0.999 * rng.random()
        assert abs(Q.renyi_entropy(math.hypot(s1, s2), s3, alpha) - _by_definition(s1, s2, s3, alpha)) < 1e-8


def test_alpha_entropy_ends():
    # alpha 2 at no bias is eq. 7; alpha 1 is the conditional von Neumann entropy; a pure state is a coin
    for c in (0.0, 0.3, 0.9, 0.995, 1.0):
        assert abs(Q.renyi_entropy(c, 0.0, 2.0) - Q.coherence_min_entropy(c)) < 1e-12
        assert abs(Q.renyi_entropy(c, 0.0, 1.0) - (1.0 - Q._h2((1.0 + c) / 2.0))) < 1e-12
        assert abs(Q.renyi_entropy(c, 0.0, 1.0 + 1e-7) - Q.renyi_entropy(c, 0.0, 1.0)) < 1e-5
    for s3 in (0.0, 0.2, 0.8):
        p0, p1 = (1 + s3) / 2, (1 - s3) / 2
        for alpha in (1.3, 2.0):
            coin = math.log2(p0 ** alpha + p1 ** alpha) / (1 - alpha)
            assert abs(Q.renyi_entropy(math.sqrt(1 - s3 * s3), s3, alpha) - coin) < 1e-7
    assert Q.renyi_entropy(0.0, 0.3, 1.5) == 0.0 and Q.renyi_entropy(0.9, 1.0, 1.5) == 0.0


def test_alpha_entropy_rises_with_coherence_and_first_rises_with_bias():
    for alpha in (1.0, 1.3, 2.0):
        by_c = [Q.renyi_entropy(c, 0.1, alpha) for c in np.linspace(0.0, 0.99, 50)]
        assert all(b >= a - 1e-12 for a, b in zip(by_c, by_c[1:]))
        assert Q.renyi_entropy(0.99, 0.05, alpha) > Q.renyi_entropy(0.99, 0.0, alpha)
        assert Q.renyi_entropy(0.99, 0.9, alpha) < Q.renyi_entropy(0.99, 0.0, alpha)


def test_worst_case_is_the_least_over_the_region():
    rng = np.random.default_rng(2)
    for _ in range(40):
        c_low, s_low = rng.random(), 0.5 * rng.random()
        s_high = min(0.999, s_low + 0.4 * rng.random())
        alpha = (1.0, 1.2, 2.0)[rng.integers(3)]
        grid = min(Q.renyi_entropy(c_low, s, alpha) for s in np.linspace(s_low, s_high, 80))
        assert Q.worst_entropy(c_low, s_low, s_high, alpha) <= grid + 1e-12
    # zero bias inside the interval is the worst case at high coherence, not the largest bias
    assert Q.worst_entropy(0.99, 0.0, 0.05, 1.0) == Q.renyi_entropy(0.99, 0.0, 1.0)


def test_shot_bound_moves_the_stated_way():
    qubit = [(0.995, 0.0, 0.03)]
    k, order = Q.shot_entropy(qubit, 20000, 2.0 ** -65)
    assert 1.0 < order <= 2.0 and 20000 * Q.coherence_min_entropy(0.995) < k < 20000 * Q.renyi_entropy(0.995, 0.0, 1.0)
    assert Q.shot_entropy(qubit, 200000, 2.0 ** -65)[0] / 200000 > k / 20000          # more shots, higher rate
    assert Q.shot_entropy(qubit, 20000, 2.0 ** -129)[0] < k                            # smaller epsilon, fewer bits
    assert Q.shot_entropy([(0.99, 0.0, 0.03)], 20000, 2.0 ** -65)[0] < k               # less coherence, fewer bits
    three = Q.shot_entropy(qubit * 3, 20000, 2.0 ** -65)[0]                             # the smoothing is paid once
    assert 3 * k < three < 60000 * Q.renyi_entropy(0.995, 0.0, 1.0)
    assert Q.shot_entropy(qubit, 50, 2.0 ** -65) == (0.0, None)                         # too few shots to certify


def test_shot_bound_is_below_the_exact_value_for_fair_coins():
    # n fair bits with nothing held by the adversary have exactly n bits of smooth min-entropy at small epsilon
    for n in (500, 5000):
        k, _ = Q.shot_entropy([(1.0, 0.0, 0.0)], n, 2.0 ** -33)
        assert n - 70 < k <= n


def _source(n, coherence=0.995, qubits=2, seed=3):
    shots, _ = Q.synthetic_qrng(n, qubits, readout_F=(0.99, 0.98), seed=seed)
    x, y = Q.synthetic_tomography(4 * n, qubits, coherence=coherence, seed=seed)
    return Q.QuantumArraySource(shots, readout=[M] * qubits, calibration_shots=4 * n, x_shots=x, y_shots=y)


def test_attestation_plans_from_the_shots_taken_together():
    att = attest(_source(20000), include_predictors=False, run_iid=True)
    plan = att.extraction_plan()
    one = es.extraction_plan(40000, att.attested_min_entropy)
    assert plan.order is not None and plan.output_bits > one.output_bits
    k, order, term = att.shot_bound(2.0 ** -65)
    assert plan.output_bits == math.floor(k - 128) and order == plan.order and term == plan.term
    # the classical values and the joint bound still cap it
    cap = min(sum(q["h_classical"] for q in att.per_qubit) / 2, att.independence["joint_rate_bound"])
    assert k <= 40000 * cap + 1e-9 and all(q["h_long_run"] > q["h_coherence"] for q in att.per_qubit)
    # without X and Y shots there is no such bound
    plain = attest(Q.QuantumArraySource(att_shots := _source(4000).shots(), readout=[M, M], calibration_shots=16000),
                   include_predictors=False)
    assert plain.shot_bound(2.0 ** -65) is None and plain.extraction_plan().order is None and att_shots.shape == (4000, 2)


def test_certificate_and_extraction_use_the_plan():
    src = _source(20000)
    cert = es.certify_source(src, select_qubits=False)
    assert "taken together, set by" in str(cert) and cert.plan.order is not None
    bits, cert2, seed = es.extract_from_source(src, seed=5)
    assert len(bits) == cert2.plan.output_bits == cert.plan.output_bits and set(np.unique(bits)) <= {0, 1}
    assert abs(float(bits.mean()) - 0.5) < 0.02


def test_saved_ibm_fez_record_gives_its_output_length():
    # job db0rsb3id5ic73eq8em0: the output length in the record, recomputed from its per-shot outcomes
    import json
    import os
    from qrng_attest.attest import pairwise_joint_bound
    base = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "hardware_coherence_ibm_fez")
    shots, meta = np.load(base + "_shots.npz"), json.load(open(base + ".json"))
    qubits = []
    for i in range(4):
        Mq, n = Q.calibrate_readout(shots["prep0"][:, i], shots["prep1"][:, i])
        r = Q.attested_coherence_entropy(shots["x"][:, i], shots["y"][:, i], Mq, calib_shots=n)
        z = shots["z"][:, i].astype(int)
        qubits.append((r["coherence_low"], Q.stokes_lower(z, Mq, calib_shots=n)[1], Q.stokes_upper(z, Mq, calib_shots=n)))
    k, order = Q.shot_entropy(qubits, 20000, 2.0 ** -65)
    cap = min(sum(row["h_classical"] for row in meta["per_qubit"]) / 4, pairwise_joint_bound(shots["z"].astype(int))[0])
    assert math.floor(min(k, 80000 * cap) - 128) == meta["extractable"] == 75020
    assert abs(order - meta["alpha"]) < 1e-12 and meta["extractable"] > math.floor(80000 * meta["attested"] - 126)


def test_saved_rounds_record_gives_its_bounds_and_output_length():
    # job db1693rid5ic73eqnak0: the five rounds pooled, each qubit's bound, the attested value and the output length
    import json
    import os
    from qrng_attest.attest import pairwise_joint_bound
    base = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results",
                        "hardware_coherence_rounds_ibm_fez")
    shots, meta = np.load(base + "_shots.npz"), json.load(open(base + ".json"))
    assert meta["plan"] == "open" and meta["pricing_type"] == "free"
    assert meta["usage"]["quantum_seconds"] <= meta["cap_s"] and meta["watch"] == "DONE"
    data = {label: np.concatenate([shots[f"{label}_{r}"] for r in range(meta["rounds"])]).astype(int)
            for label in ("prep0", "prep1", "z", "x", "y")}
    assert data["z"].shape == (20000, 12) and data["x"].shape == (100000, 12)
    qubits, bounds = [], []
    for i, row in enumerate(meta["per_qubit"]):
        Mq, n = Q.calibrate_readout(data["prep0"][:, i], data["prep1"][:, i])
        r = Q.attested_coherence_entropy(data["x"][:, i], data["y"][:, i], Mq, calib_shots=n)
        assert abs(r["h_coherence"] - row["h_coherence"]) < 1e-12 and 0.99 < r["coherence_low"] <= r["coherence"]
        z = data["z"][:, i]
        qubits.append((r["coherence_low"], Q.stokes_lower(z, Mq, calib_shots=n)[1], Q.stokes_upper(z, Mq, calib_shots=n)))
        bounds.append(r["h_coherence"])
    assert abs(min(bounds) - meta["attested"]) < 1e-9 and abs(meta["attested"] - 0.8149) < 5e-5
    k, order = Q.shot_entropy(qubits, 20000, 2.0 ** -65)
    classical = 20000 * sum(row["h_classical"] for row in meta["per_qubit"])
    joint = 240000 * pairwise_joint_bound(data["z"])[0]
    assert classical < k < joint and abs(order - meta["alpha"]) < 1e-12
    assert math.floor(classical - 128) == meta["extractable"] == 229744
    # no qubit's coherence moves between rounds by more than three parts in a thousand
    assert max(max(c) - min(c) for c in meta["coherence_by_round"]) < 0.003

