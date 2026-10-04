# the coherence bound: eq. 7 of Fiorentino et al., its estimate from X and Y basis shots, and the attestation
import contextlib
import io
import math

import numpy as np

import qrng_attest as es
from _harness import skip
from qrng_attest import quantum as Q
from qrng_attest.attest import attest, qubit_subset

M = Q.assignment_matrix(0.99, 0.98)


def _z_shots(n, seed, qubits=1):
    # Z outcomes of a state with S3 = 0 through the readout M
    rng = np.random.default_rng([seed, 2])
    true = rng.integers(0, 2, (n, qubits))
    return np.where(rng.random((n, qubits)) < np.where(true == 0, 0.01, 0.02), 1 - true, true)


def test_coherence_min_entropy_values():
    assert Q.coherence_min_entropy(1.0) == 1.0 and Q.coherence_min_entropy(0.0) == 0.0
    assert abs(Q.coherence_min_entropy(0.99) - 0.8096) < 5e-5 and abs(Q.coherence_min_entropy(0.9) - 0.4781) < 5e-5
    grid = [Q.coherence_min_entropy(c) for c in np.linspace(0, 1, 101)]
    assert all(b > a for a, b in zip(grid, grid[1:]))


def test_eq7_is_the_guess_of_an_adversary_holding_a_purification():
    # -log2 of the Helstrom guessing probability between the adversary's two states after the Z outcome
    rng = np.random.default_rng(35)
    for _ in range(200):
        v = rng.normal(size=2) + 1j * rng.normal(size=2)
        v /= np.linalg.norm(v)
        t = rng.random()
        rho = t * np.outer(v, v.conj()) + (1 - t) * np.eye(2) / 2
        lam, vec = np.linalg.eigh(rho)
        lam = np.clip(lam, 0.0, None)
        tau = [np.outer(np.sqrt(lam) * vec[z, :], (np.sqrt(lam) * vec[z, :]).conj()) for z in (0, 1)]
        guess = 0.5 * (1.0 + float(np.abs(np.linalg.eigvalsh(tau[0] - tau[1])).sum()))
        assert abs(Q.coherence_min_entropy(2.0 * abs(rho[0, 1])) + math.log2(guess)) < 1e-9


def test_population_value_is_at_or_above_eq7():
    # |S3| <= sqrt(1 - c^2) for a physical state, so -log2 max(P0, P1) cannot be below eq. 7
    rng = np.random.default_rng(1)
    for _ in range(500):
        s = rng.normal(size=3)
        s *= rng.random() / np.linalg.norm(s)
        assert Q.geometric_min_entropy((1 + s[2]) / 2) >= Q.coherence_min_entropy(math.hypot(s[0], s[1])) - 1e-12


def test_the_bound_is_at_or_below_the_truth():
    for c in (1.0, 0.99, 0.9, 0.5, 0.0):
        truth = Q.coherence_min_entropy(c)
        for seed in range(20):
            x, y = Q.synthetic_tomography(20000, 1, coherence=c, seed=seed)
            r = Q.attested_coherence_entropy(x[:, 0], y[:, 0], M, calib_shots=20000)
            assert r["h_coherence"] <= truth + 1e-12, (c, seed, r)
            assert r["coherence_low"] <= r["coherence"] + 1e-12


def test_a_classical_coin_is_attested_at_zero():
    # uniform Z outcomes with no coherence: the X and Y shots are uniform too
    shots = _z_shots(20000, 0)
    x, y = Q.synthetic_tomography(20000, 1, coherence=0.0, seed=0)
    with_basis = attest(Q.QuantumArraySource(shots, readout=[M], calibration_shots=20000, x_shots=x, y_shots=y),
                        include_predictors=False)
    without = attest(Q.QuantumArraySource(shots, readout=[M], calibration_shots=20000), include_predictors=False)
    assert with_basis.attested_min_entropy < 0.01 and not with_basis.passed
    assert with_basis.per_qubit[0]["h_quantum"] < 0.01 < 0.9 < without.per_qubit[0]["h_quantum"]
    assert any("assumes a pure source state" in f for f in without.describe["flags"])
    assert not any("assumes a pure source state" in f for f in with_basis.describe["flags"])


def test_all_zero_x_shots_do_not_certify_full_coherence():
    # 200 shots, all 0: the lower bound on P(0) is 200 / (200 + z^2), not 1
    r = Q.attested_coherence_entropy(np.zeros(200, int))
    z2 = Q.Z_CONF ** 2
    assert abs(r["coherence_low"] - (2 * 200 / (200 + z2) - 1)) < 1e-12 and r["h_coherence"] < 0.8
    assert Q.attested_coherence_entropy(np.zeros(200000, int))["h_coherence"] > 0.98


def test_each_input_moves_the_bound_the_stated_way():
    x, y = Q.synthetic_tomography(40000, 1, coherence=0.995, seed=1)
    base = Q.attested_coherence_entropy(x[:, 0], y[:, 0], M, calib_shots=40000)["h_coherence"]
    assert Q.attested_coherence_entropy(x[:2000, 0], y[:, 0], M, calib_shots=40000)["h_coherence"] < base
    assert Q.attested_coherence_entropy(x[:, 0], y[:, 0], M, calib_shots=500)["h_coherence"] < base
    assert Q.attested_coherence_entropy(x[:, 0], None, M, calib_shots=40000)["h_coherence"] <= base
    # coherence along Y is credited only with the Y shots
    rng = np.random.default_rng(7)
    flat, along_y = rng.integers(0, 2, 40000), (rng.random(40000) > 0.995).astype(int)
    assert Q.attested_coherence_entropy(flat, None)["h_coherence"] < 0.01
    assert Q.attested_coherence_entropy(flat, along_y)["h_coherence"] > 0.7


def test_reset_error_is_not_credited_as_coherence():
    # a qubit that starts in 1 with probability p, read perfectly: the calibration reads F0 = 1 - p, the state has 1 - 2p
    rng = np.random.default_rng(4)
    n = 50000
    for p in (0.002, 0.01, 0.05):
        start = rng.random((3, n)) < p
        Mq, nc = Q.calibrate_readout(start[0].astype(int), (~start[1]).astype(int))
        r = Q.attested_coherence_entropy(start[2].astype(int), rng.integers(0, 2, n), Mq, calib_shots=nc)
        assert r["coherence_low"] <= 1 - 2 * p and r["h_coherence"] <= Q.coherence_min_entropy(1 - 2 * p), (p, r)


def test_attestation_takes_the_smaller_of_the_two_values():
    shots = _z_shots(20000, 3, qubits=2)
    x, y = Q.synthetic_tomography(20000, 2, coherence=0.9, seed=3)
    att = attest(Q.QuantumArraySource(shots, readout=[M, M], calibration_shots=20000, x_shots=x, y_shots=y),
                 include_predictors=False)
    for q in att.per_qubit:
        assert q["h_quantum"] == min(q["h_population"], q["h_coherence"]) and q["h_coherence"] < 0.48
    assert att.attested_min_entropy < 0.48
    assert any("coherence bound sets the attested rate" in a for a in es.advice(att))
    assert es.remediation(att) == ["no action required; re-attest after the next calibration cycle"] and att.passed


def test_subset_search_carries_the_basis_shots():
    shots = _z_shots(8000, 4, qubits=3)
    x, y = Q.synthetic_tomography(8000, 3, coherence=0.99, seed=4)
    x[:, 1] = np.random.default_rng(4).integers(0, 2, 8000)             # qubit 1 has lost its coherence
    y[:, 1] = np.random.default_rng(5).integers(0, 2, 8000)
    sub = qubit_subset(Q.QuantumArraySource(shots, readout=[M] * 3, calibration_shots=8000, x_shots=x, y_shots=y))
    assert 1 in sub["dropped"] and sub["best_rate"] > 0.5 > sub["base_rate"]


def test_demo_reports_the_coherence():
    from qrng_attest.__main__ import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        main(["demo", "--shots", "4000", "--qubits", "2"])
        main(["demo", "--shots", "4000", "--qubits", "2", "--no-tomography"])
    text = out.getvalue()
    assert "source coherence at least" in text and "assumes a pure source state" in text


def test_x_basis_circuit_keeps_its_preparation():
    # without the barrier the transpiler cancels H H and the circuit would measure |0> for any preparation
    try:
        from qiskit import transpile
        from qiskit.quantum_info import Statevector
    except ImportError:
        skip("requires qiskit")
    from qrng_attest.qiskit_bridge import source_circuit
    for basis, want in (("z", 0.0), ("x", 1.0), ("y", 0.0)):
        qc = source_circuit(2, basis)
        p0 = Statevector(qc.remove_final_measurements(inplace=False)).probabilities([0])[0]
        assert abs(2 * p0 - 1 - want) < 1e-12, basis
    t = transpile(source_circuit(1, "x"), basis_gates=["sx", "rz", "x"], optimization_level=3)
    assert sum(1 for i in t.data if i.operation.name in ("sx", "rz", "x")) > 0


def test_backend_simulation_carries_x_and_y_shots():
    try:
        import qiskit_aer  # noqa: F401
        from qiskit_ibm_runtime.fake_provider import FakeManilaV2
    except ImportError:
        skip("requires qiskit-aer and qiskit-ibm-runtime")
    from qrng_attest.qiskit_bridge import quantum_source_from_backend
    src = quantum_source_from_backend(FakeManilaV2(), n_shots=4000, qubits=[0, 1], seed=5)
    basis = src.basis_shots()
    assert basis["x"].shape == basis["y"].shape == (4000, 2)
    att = attest(src, include_predictors=False)
    for q in att.per_qubit:
        assert 0.9 < q["coherence"] <= 1.0 and q["h_coherence"] <= q["h_population"] + 1e-12
    assert quantum_source_from_backend(FakeManilaV2(), n_shots=500, qubits=[0], tomography=False).basis_shots() is None


def test_backend_command_finds_a_fake_backend_by_name():
    # "manila" is also the name of a module of the fake provider, which is not a backend
    try:
        import qiskit_ibm_runtime.fake_provider  # noqa: F401
    except ImportError:
        skip("requires qiskit-ibm-runtime")
    from qrng_attest.__main__ import _load_backend
    assert _load_backend("manila").num_qubits == 5 and _load_backend("FakeManilaV2").num_qubits == 5
    assert _load_backend("no_such_device") is None


def test_hardware_script_plans_on_a_fake_backend():
    # uncoupled qubits, every circuit read on them, and the basis changes still in the transpiled circuits
    try:
        import qiskit_ibm_runtime.fake_provider as fp
        from qiskit.transpiler import generate_preset_pass_manager
    except ImportError:
        skip("requires qiskit-ibm-runtime")
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
    import coherence_hardware as ch
    backend = fp.FakeManilaV2()
    qubits = ch.pick_qubits(backend, 2)
    edges = {tuple(sorted(e)) for e in backend.coupling_map.get_edges()}
    assert len(qubits) == 2 and tuple(sorted(qubits)) not in edges
    pm = generate_preset_pass_manager(optimization_level=1, backend=backend, initial_layout=qubits,
                                      scheduling_method="alap")
    tqcs = [pm.run(c) for c in ch.circuits(2)]
    assert all(ch.measured(t) == qubits for t in tqcs)
    count = dict(zip(ch.LABELS, (ch.pulses(t) for t in tqcs)))
    assert count["prep0"] == 0 and count["x"] > count["z"] and count["y"] > count["z"]


def test_rounds_script_pools_and_tabulates():
    # two rounds of a synthetic source: the pooled certificate, the coherence of a round, the bound after each round
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
    import coherence_hardware as ch
    import coherence_rounds as cr
    rng = np.random.default_rng(3)
    got = {}
    for r in range(2):
        got[f"prep0_{r}"] = (rng.random((4000, 2)) > 0.99).astype(int)
        got[f"prep1_{r}"] = (rng.random((4000, 2)) < 0.98).astype(int)
        got[f"z_{r}"] = _z_shots(1500, 30 + r, 2)
        got[f"x_{r}"], got[f"y_{r}"] = Q.synthetic_tomography(4000, 2, coherence=0.99, seed=20 + r)
    text, data = cr.analyse(got, [5, 7], 2, "synthetic", "synthetic")
    cert, _ = ch.certify(cr.pooled(got, 2), "synthetic", "synthetic")
    assert abs(cert.rate - data["attested"]) < 1e-12 and "coherence measured in each round" in text
    for i in range(2):
        Mq, n = Q.calibrate_readout(got["prep0_0"][:, i], got["prep1_0"][:, i])
        first = Q.attested_coherence_entropy(got["x_0"][:, i], got["y_0"][:, i], Mq, calib_shots=n)
        assert abs(data["coherence_by_round"][i][0] - first["coherence"]) < 1e-12
        assert abs(data["bound_after_rounds"][i][0] - first["h_coherence"]) < 1e-12
        assert abs(data["bound_after_rounds"][i][1] - data["per_qubit"][i]["h_coherence"]) < 1e-12
    assert [n for _label, _r, n in cr.pubs()] == [20000, 20000, 4000, 20000, 20000] * 5


def test_rounds_script_stops_each_unsafe_submission():
    # one bad plan per safety, each stopped for its own reason and no other; the plan as made is not stopped
    import os
    import sys
    import time
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
    import coherence_rounds as cr
    now = 1_800_000_000
    stamp = lambda t: time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t))
    good = {"made_utc": stamp(now - 300), "pubs": [list(x) for x in cr.pubs()], "qubits": list(range(cr.QUBITS)),
            "estimate_s": 106.0, "cap_s": 169}
    usage = {"usage_remaining_seconds": 478, "usage_limit_reached": False}
    assert cr.safeties(good, usage, now) == []
    cases = {"old": ({**good, "made_utc": stamp(now - 7200)}, usage, "minutes old"),
             "larger": ({**good, "pubs": good["pubs"] + [["x", 5, 20000]]}, usage, "not this script's job"),
             "other shots": ({**good, "pubs": [[l, r, 2 * n] for l, r, n in good["pubs"]]}, usage, "not this script's job"),
             "fewer qubits": ({**good, "qubits": [0, 1]}, usage, "not this script's job"),
             "estimate": ({**good, "estimate_s": 150.0, "cap_s": 170}, usage, "the estimate"),
             "cap": ({**good, "cap_s": 400}, {**usage, "usage_remaining_seconds": 600}, "the cap"),
             "reserve": (good, {**usage, "usage_remaining_seconds": 250}, "would leave less"),
             "limit reached": (good, {**usage, "usage_limit_reached": True}, "would leave less")}
    for name, (plan, use, phrase) in cases.items():
        reasons = cr.safeties(plan, use, now)
        assert len(reasons) == 1 and phrase in reasons[0], (name, reasons)
    assert not cr.overrun(None, now, 169) and not cr.overrun(now - 180, now, 169) and cr.overrun(now - 190, now, 169)
    assert cr._seconds("2026-10-04T05:12:33.123456Z") == cr._seconds("2026-10-04T05:12:33Z")


def test_saved_ibm_fez_record_gives_its_coherence_bounds():
    # job db0rsb3id5ic73eq8em0: the per-qubit bounds in the record, recomputed from its per-shot outcomes
    import json
    import os
    base = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results",
                        "hardware_coherence_ibm_fez")
    shots, meta = np.load(base + "_shots.npz"), json.load(open(base + ".json"))
    assert meta["plan"] == "open" and meta["pricing_type"] == "free" and shots["z"].shape == (20000, 4)
    bounds = []
    for i, row in enumerate(meta["per_qubit"]):
        Mq, n = Q.calibrate_readout(shots["prep0"][:, i], shots["prep1"][:, i])
        r = Q.attested_coherence_entropy(shots["x"][:, i], shots["y"][:, i], Mq, calib_shots=n)
        assert abs(r["coherence_low"] - row["coherence_low"]) < 1e-12 and abs(r["h_coherence"] - row["h_coherence"]) < 1e-12
        assert 0.98 < r["coherence_low"] <= r["coherence"] <= 1.0
        bounds.append(r["h_coherence"])
    assert abs(min(bounds) - meta["attested"]) < 1e-9 and abs(meta["attested"] - 0.7653) < 5e-5

