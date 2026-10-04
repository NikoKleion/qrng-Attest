# what a first reader meets: the shapes a source takes, the verdict, the command line, and the examples in USAGE.md
import contextlib
import io
import os
import re
import sys
import tempfile

import numpy as np

import qrng_attest as qa
from _harness import skip
from qrng_attest import quantum as Q
from qrng_attest.__main__ import _load_backend, main

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M = Q.assignment_matrix(0.99, 0.98)


def _run(argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = main(argv)
    return code, out.getvalue()


def test_one_readout_matrix_serves_every_qubit_and_wrong_shapes_say_so():
    z, _ = Q.synthetic_qrng(3000, 2, seed=1)
    x, y = Q.synthetic_tomography(6000, 2, seed=1)
    one = qa.attest(qa.QuantumArraySource(z, readout=M, calibration_shots=6000, x_shots=x, y_shots=y), include_predictors=False)
    each = qa.attest(qa.QuantumArraySource(z, readout=[M, M], calibration_shots=6000, x_shots=x, y_shots=y),
                     include_predictors=False)
    assert one.attested_min_entropy == each.attested_min_entropy
    for bad, phrase in ((dict(readout=[M]), "1 matrices for 2 qubits"),
                        (dict(readout=M, x_shots=x[:, :1], y_shots=y), "x_shots has 1 columns for 2 qubits")):
        try:
            qa.QuantumArraySource(z, **bad)
        except ValueError as e:
            assert phrase in str(e), e
        else:
            raise AssertionError(bad.keys())


def test_a_source_without_x_and_y_shots_is_not_passed():
    # the Z shots of a source with no coherence look like a fair coin; without X and Y shots nothing is attested
    z, _ = Q.synthetic_qrng(6000, 2, seed=2)
    src = qa.QuantumArraySource(z, readout=M, calibration_shots=6000)
    cert = qa.certify_source(src, select_qubits=False)
    text = str(cert)
    assert not cert.passed and "FLAGGED" in text and "not extracted from while FLAGGED" in text
    assert any("no X or Y basis shots" in a for a in cert.actions)
    try:
        qa.extract_from_source(src)
    except RuntimeError:
        pass
    else:
        raise AssertionError("extracted from a source that is not attested")
    x, y = Q.synthetic_tomography(6000, 2, coherence=0.0, seed=2)
    none = qa.certify_source(qa.QuantumArraySource(z, readout=M, calibration_shots=6000, x_shots=x, y_shots=y),
                             select_qubits=False)
    assert none.rate == 0.0 and not none.passed and any("no coherence is certified" in a for a in none.actions)
    assert "no number of raw bits" in str(none)


def test_failed_checks_are_actions_and_the_rest_is_advice():
    z, _ = Q.synthetic_qrng(6000, 2, seed=3)
    x, y = Q.synthetic_tomography(24000, 2, coherence=0.99, seed=3)
    cert = qa.certify_source(qa.QuantumArraySource(z, readout=M, calibration_shots=24000, x_shots=x, y_shots=y),
                             select_qubits=False)
    assert cert.passed and cert.actions == ["no action required; re-attest after the next calibration cycle"]
    assert cert.advice and all(line in str(cert) for line in cert.advice) and "  advice: " in str(cert)
    bad, _ = Q.synthetic_qrng(6000, 2, crosstalk=0.6, seed=3)
    flagged = qa.certify_source(qa.QuantumArraySource(bad, readout=M, calibration_shots=24000, x_shots=x, y_shots=y),
                                select_qubits=False)
    assert not flagged.passed and any("correlated" in a for a in flagged.actions)
    assert not any("correlated" in a for a in flagged.advice)


def test_capture_command_takes_its_documented_flags():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "qrng.bin")
        np.random.default_rng(4).integers(0, 256, 1500, dtype=np.uint8).tofile(path)
        for flags in ([], ["--no-iid"], ["--iid-full"], ["--bytes", "--no-iid"]):
            code, out = _run(["capture", path] + flags)
            assert code == 0 and "min-entropy =" in out, flags
        code, out = _run(["capture", os.path.join(tmp, "absent.bin")])
        assert code == 1 and "no such file" in out


def test_a_byte_file_names_the_estimator_that_sets_its_value():
    rng = np.random.default_rng(5)
    S = (rng.random(6000) < 0.03) * rng.integers(1, 256, 6000)
    a = qa.assess(S.astype(int), run_iid=False)
    word = a.meta["word_size"]
    track, _name = a.driver.split(":")
    value = a.scores[a.driver] * (word if track == "bitstring" else 1.0)
    assert abs(value - a.min_entropy) < 1e-9, (a.driver, value, a.min_entropy)


def test_backend_names_resolve_to_device_models():
    try:
        import qiskit_ibm_runtime.fake_provider  # noqa: F401
    except ImportError:
        skip("requires qiskit-ibm-runtime")
    assert _load_backend("manila").num_qubits == 5
    assert _load_backend("fez").num_qubits == _load_backend("ibm_fez").num_qubits == 156
    assert _load_backend("nosuchdevice") is None
    code, out = _run(["backend", "nosuchdevice"])
    assert code == 1 and "fake backends" in out and "certify_backend(backend)" in out


def test_rounds_script_keeps_a_saved_record():
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import time
    import coherence_rounds as cr
    now = 1_800_000_000
    plan = {"made_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(now - 60)), "pubs": [list(x) for x in cr.pubs()],
            "qubits": list(range(cr.QUBITS)), "estimate_s": 106.0, "cap_s": 169}
    usage = {"usage_remaining_seconds": 478, "usage_limit_reached": False}
    saved = os.path.join(ROOT, "results", "hardware_coherence_rounds_ibm_fez_job.json")
    assert os.path.exists(saved)
    reasons = cr.safeties(plan, usage, now, record=saved)
    assert len(reasons) == 1 and "exists" in reasons[0]
    assert cr.safeties(plan, usage, now, record=saved + ".absent") == []


def _blocks(heading):
    text = open(os.path.join(ROOT, "USAGE.md"), encoding="utf-8").read()
    section = text[text.index(heading):]
    section = section[:section.index("\n## ", 5)]
    return re.findall(r"```python\n(.*?)```", section, re.S)


def test_the_arrays_example_of_the_usage_guide_runs_as_written():
    (block,) = _blocks("## From arrays of shots to uniform bits")
    scope = {}
    with contextlib.redirect_stdout(io.StringIO()) as out:
        exec(block, scope)
    assert "QRNG certificate" in out.getvalue() and scope["cert"].passed
    assert len(scope["bits"]) == scope["cert"].plan.output_bits > 15000


def test_the_self_testing_example_of_the_usage_guide_runs_as_written():
    try:
        import qiskit_aer  # noqa: F401
        import qiskit_ibm_runtime.fake_provider  # noqa: F401
    except ImportError:
        skip("requires qiskit-aer and qiskit-ibm-runtime")
    (block,) = _blocks("## Self-testing bound")
    scope = {}
    with contextlib.redirect_stdout(io.StringIO()) as out:
        exec(block, scope)
    printed = out.getvalue().split()
    assert len(scope["circuits"]) == 8 and sorted(scope["circuits"])[0] == (0, 0)
    assert 0.0 < float(printed[2]) < 0.2285 and 0.2 < float(printed[3]) < 0.2285
