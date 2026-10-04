# the source in the Z, X and Y bases with its readout calibration, on an IBM backend, in one job on one set of qubits
#   python scripts/coherence_hardware.py dry [backend]   plan, estimate and rehearsal; nothing is submitted
#   python scripts/coherence_hardware.py submit          checks the plan and the allowance, then submits
#   python scripts/coherence_hardware.py collect         reads the finished job and certifies
import json
import math
import os
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
OUT = os.path.join(REPO, "build", "coherence_dry_run")
QUBITS = 4
SHOTS = 20000
CAP_SECONDS = 90           # hard cap on the job's execution time
REHEARSAL_SHOTS = 2000     # per circuit, under the backend's noise model, at the dry run
LABELS = ("prep0", "prep1", "z", "x", "y")


def guards(svc, need_seconds):
    insts = svc.instances()
    if len(insts) != 1:
        raise SystemExit(f"expected exactly one instance, found {len(insts)}; stopping")
    inst = insts[0]
    if inst.get("plan") != "open" or inst.get("pricing_type") != "free":
        raise SystemExit(f"instance is not the free open plan ({inst.get('plan')}, {inst.get('pricing_type')}); stopping")
    u = svc.usage()
    if u.get("usage_limit_reached") or u.get("usage_remaining_seconds", 0) < need_seconds:
        raise SystemExit(f"allowance too low: {u.get('usage_remaining_seconds')} s left, {need_seconds} s needed; stopping")
    return inst, u


def pick_qubits(backend, n):
    # the n qubits of least readout plus sx error, no two of them coupled
    target = backend.target
    neighbours = {q: set() for q in range(target.num_qubits)}
    for a, b in backend.coupling_map.get_edges():
        neighbours[a].add(b)
        neighbours[b].add(a)
    ranked = []
    for q in range(target.num_qubits):
        m, sx, props = target["measure"].get((q,)), target["sx"].get((q,)), target.qubit_properties[q]
        if m is None or sx is None or props is None or m.error is None or sx.error is None or not props.t2:
            continue
        ranked.append((m.error + sx.error, q))
    chosen = []
    for _err, q in sorted(ranked):
        if all(q not in neighbours[p] for p in chosen):
            chosen.append(q)
        if len(chosen) == n:
            return chosen
    raise SystemExit(f"fewer than {n} usable uncoupled qubits; stopping")


def circuits(nq):
    from qiskit import QuantumCircuit
    from qrng_attest.qiskit_bridge import source_circuit
    c0 = QuantumCircuit(nq, nq)
    c0.measure(range(nq), range(nq))
    c1 = QuantumCircuit(nq, nq)
    for i in range(nq):
        c1.x(i)
    c1.measure(range(nq), range(nq))
    return [c0, c1] + [source_circuit(nq, basis) for basis in ("z", "x", "y")]


def measured(tqc):
    # physical qubit read into each classical bit
    idx = lambda q: tqc.find_bit(q).index
    return [p for _c, p in sorted((tqc.find_bit(i.clbits[0]).index, idx(i.qubits[0])) for i in tqc.data
                                  if i.operation.name == "measure")]


def pulses(tqc):
    return sum(1 for i in tqc.data if i.operation.name in ("sx", "x"))


def calibration(target, qubits):
    out = {}
    for q in qubits:
        p = target.qubit_properties[q]
        out[str(q)] = {"t1": p.t1, "t2": p.t2, "readout_error": target["measure"][(q,)].error,
                       "measure_duration": target["measure"][(q,)].duration, "sx_error": target["sx"][(q,)].error}
    return out


def dry(name=None):
    from qiskit import qpy
    from qiskit.transpiler import generate_preset_pass_manager
    from qiskit_ibm_runtime import QiskitRuntimeService

    svc = QiskitRuntimeService()
    inst, usage = guards(svc, 1)
    print(f"instance {inst.get('name')}: plan {inst.get('plan')}, pricing {inst.get('pricing_type')}; "
          f"{usage['usage_remaining_seconds']} s of {usage['usage_limit_seconds']} left this period")
    backend = svc.backend(name) if name else svc.least_busy(operational=True, simulator=False)
    status, target = backend.status(), backend.target
    print(f"{backend.name}: operational {status.operational}, pending jobs {status.pending_jobs}, "
          f"{target.num_qubits} qubits")

    qubits = pick_qubits(backend, QUBITS)
    cal = calibration(target, qubits)
    print(f"qubits {qubits}")
    for q in qubits:
        c = cal[str(q)]
        print(f"  {q}: readout error {c['readout_error']:.4f}, sx error {c['sx_error']:.1e}, "
              f"T1 {c['t1'] * 1e6:.0f} us, T2 {c['t2'] * 1e6:.0f} us")

    pm = generate_preset_pass_manager(optimization_level=1, backend=backend, initial_layout=qubits,
                                      scheduling_method="alap")
    tqcs = [pm.run(c) for c in circuits(QUBITS)]
    for label, t in zip(LABELS, tqcs):
        assert measured(t) == qubits, f"circuit {label} reads {measured(t)}, not {qubits}"
    count = dict(zip(LABELS, (pulses(t) for t in tqcs)))
    assert count["x"] > count["z"] and count["y"] > count["z"], f"a basis change was cancelled: {count}"
    print(f"all {len(tqcs)} circuits on those qubits; sx and x pulses per circuit {count}")

    durs = [t.estimate_duration(target, unit="s") for t in tqcs]
    rep = backend.configuration().default_rep_delay
    rep = rep * 1e-6 if rep and rep > 1e-3 else (rep or 250e-6)
    est = sum(SHOTS * (d + rep) for d in durs)
    need = math.ceil(2 * est)
    print(f"scheduled durations {[round(d * 1e6, 2) for d in durs]} us; rep delay {rep * 1e6:.0f} us; "
          f"{SHOTS} shots per circuit")
    print(f"estimate: {est:.0f} s of shots plus overhead; submit needs {need} s left and caps the job at "
          f"{CAP_SECONDS} s")
    guards(svc, need)

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "circuits.qpy"), "wb") as fh:
        qpy.dump(tqcs, fh)
    plan = {"backend": backend.name, "qubits": qubits, "shots": SHOTS, "labels": list(LABELS),
            "durations_s": durs, "rep_delay_s": rep, "estimate_s": est, "need_s": need, "cap_s": CAP_SECONDS,
            "calibration": cal, "made_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(plan, open(os.path.join(OUT, "plan.json"), "w"), indent=1)

    try:
        from qiskit_aer import AerSimulator
    except ImportError:
        print("qiskit-aer is not installed; no rehearsal")
    else:
        print()
        print(f"rehearsal of the same circuits under the Aer noise model of {backend.name}, {REHEARSAL_SHOTS} shots each")
        sim = AerSimulator.from_backend(backend)
        got = {}
        for k, (label, t) in enumerate(zip(LABELS, tqcs)):
            mem = sim.run(t, shots=REHEARSAL_SHOTS, memory=True,
                          seed_simulator=7 + k * (REHEARSAL_SHOTS + 1)).result().get_memory()
            got[label] = np.array([[int(x) for x in s.replace(" ", "")[::-1]] for s in mem], int)
        print(certify(got, f"{backend.name} model", "synthetic")[0])
    print()
    print(f"saved to {OUT}; nothing was submitted")


def certify(got, name, provenance):
    from qrng_attest.certify import certify_source
    from qrng_attest.quantum import QuantumArraySource, calibrate_readout
    Ms, n_cal = [], len(got["prep0"])
    for i in range(got["z"].shape[1]):
        M, n = calibrate_readout(got["prep0"][:, i], got["prep1"][:, i])
        Ms.append(M)
        n_cal = min(n_cal, n)
    src = QuantumArraySource(got["z"], readout=Ms, calibration_shots=n_cal, name=f"{name} QRNG",
                             provenance=provenance, x_shots=got["x"], y_shots=got["y"])
    return certify_source(src), Ms


def submit():
    from qiskit import qpy
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    import qiskit
    import qiskit_ibm_runtime

    plan = json.load(open(os.path.join(OUT, "plan.json")))
    with open(os.path.join(OUT, "circuits.qpy"), "rb") as fh:
        tqcs = list(qpy.load(fh))
    need = math.ceil(2 * plan["estimate_s"])
    base = os.path.join(REPO, "results", f"hardware_coherence_{plan['backend']}")
    if os.path.exists(base + "_job.json"):
        raise SystemExit(f"{os.path.basename(base)}_job.json exists; move the record of the earlier job before another run")
    svc = QiskitRuntimeService()
    inst, usage = guards(svc, need)
    backend = svc.backend(plan["backend"])
    if not backend.status().operational:
        raise SystemExit(f"{plan['backend']} is not operational; stopping")
    target = backend.target
    for q in plan["qubits"]:
        if target.qubit_properties[q] is None or (q,) not in target["measure"]:
            raise SystemExit(f"qubit {q} is no longer usable on {plan['backend']}; stopping")
    sampler = SamplerV2(mode=backend)
    sampler.options.default_shots = plan["shots"]
    sampler.options.max_execution_time = plan["cap_s"]
    sampler.options.dynamical_decoupling.enable = False
    sampler.options.twirling.enable_gates = False
    sampler.options.twirling.enable_measure = False
    cal = calibration(target, plan["qubits"])
    inst, usage = guards(svc, need)          # once more, immediately before the job
    job = sampler.run(tqcs)
    stub = {"job_id": job.job_id(), "backend": plan["backend"], "qubits": plan["qubits"], "shots": plan["shots"],
            "labels": plan["labels"], "submitted_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "plan": inst.get("plan"), "pricing_type": inst.get("pricing_type"),
            "estimate_s": plan["estimate_s"],
            "cap_s": plan["cap_s"], "calibration_at_submission": cal,
            "qiskit_version": qiskit.__version__, "qiskit_ibm_runtime_version": qiskit_ibm_runtime.__version__}
    base = os.path.join(REPO, "results", f"hardware_coherence_{plan['backend']}")
    json.dump(stub, open(base + "_job.json", "w"), indent=1)
    with open(base + "_circuits.qpy", "wb") as fh:
        qpy.dump(tqcs, fh)
    print(f"submitted job {stub['job_id']} to {plan['backend']}; run collect when it has finished")


def collect():
    from qiskit_ibm_runtime import QiskitRuntimeService

    plan = json.load(open(os.path.join(OUT, "plan.json")))
    base = os.path.join(REPO, "results", f"hardware_coherence_{plan['backend']}")
    stub = json.load(open(base + "_job.json"))
    if not os.path.exists(base + "_shots.npz"):
        job = QiskitRuntimeService().job(stub["job_id"])
        status = str(job.status())
        if "DONE" not in status.upper():
            print(f"job {stub['job_id']}: {status}; nothing collected")
            return
        got = {}
        for label, pub in zip(stub["labels"], job.result()):
            strings = next(iter(pub.data.values())).get_bitstrings()
            got[label] = np.array([[int(x) for x in s.replace(" ", "")[::-1]] for s in strings], np.uint8)
        np.savez_compressed(base + "_shots.npz", **got)
        stub["usage"] = job.metrics().get("usage", {})
        stub["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        json.dump(stub, open(base + "_job.json", "w"), indent=1)
    # the certificate is computed from the saved shots, so it can be recomputed without the account
    got = {k: v.astype(int) for k, v in np.load(base + "_shots.npz").items()}
    usage = stub.get("usage", {})
    cert, Ms = certify(got, stub["backend"], "measured")
    att = cert.attestation
    lines = [f"job {stub['job_id']} on {stub['backend']}, qubits {stub['qubits']}, {stub['shots']} shots per circuit, "
             f"plan {stub['plan']} and {stub['pricing_type']}, usage {usage}",
             "", f"  {'qubit':>5} {'F0':>7} {'F1':>7} {'P(0) Z':>7} {'coherence':>10} {'at least':>9} "
                 f"{'population':>11} {'bound':>7}"]
    for q, M, row, z in zip(stub["qubits"], Ms, att.per_qubit, got["z"].T):
        lines.append(f"  {q:5d} {M[0, 0]:7.4f} {M[1, 1]:7.4f} {float(np.mean(z == 0)):7.4f} {row['coherence']:10.4f} "
                     f"{row['coherence_low']:9.4f} {row['h_population']:11.4f} {row['h_coherence']:7.4f}")
    lines += ["", str(cert)]
    text = "\n".join(lines) + "\n"
    with open(base + ".txt", "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    json.dump({**stub, "attested": cert.rate, "passed": bool(cert.passed),
               "extractable": cert.plan.output_bits, "alpha": cert.plan.order,
               "per_qubit": [{k: row[k] for k in ("coherence", "coherence_low", "h_population", "h_coherence",
                                                 "h_classical", "classical_track", "h_quantum", "s3_low",
                                                 "s3_high", "h_long_run")}
                             for row in att.per_qubit],
               "readout": [[float(M[0, 0]), float(M[1, 1])] for M in Ms]},
              open(base + ".json", "w"), indent=1)
    print(text, end="")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "dry":
        dry(sys.argv[2] if len(sys.argv) > 2 else None)
    elif mode in ("submit", "collect"):
        {"submit": submit, "collect": collect}[mode]()
    else:
        print("usage: python scripts/coherence_hardware.py dry [backend] | submit | collect")
