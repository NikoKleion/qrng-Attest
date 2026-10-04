# the source in the Z, X and Y bases with its readout calibration, in rounds, on an IBM backend, in one job
#   python scripts/coherence_rounds.py dry [backend] [norehearse]   plan, estimate and rehearsal; nothing is submitted
#   python scripts/coherence_rounds.py submit                       prints the plan and the reasons to stop, if any
#   python scripts/coherence_rounds.py submit go                    submits when there is none, then watches the job
#   python scripts/coherence_rounds.py collect                      reads the finished job, certifies, tabulates by round
import calendar
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import coherence_hardware as H

OUT = os.path.join(REPO, "build", "coherence_rounds")
QUBITS = 12
ROUNDS = 5
SHOTS = {"prep0": 20000, "prep1": 20000, "z": 4000, "x": 20000, "y": 20000}      # per round
CAP_LIMIT = 180            # no plan caps the job above this many seconds
ESTIMATE_LIMIT = 117       # seconds; a plan estimated above this is not submitted
RESERVE = 120              # seconds that stay on the allowance if the job runs to its cap
PLAN_MAX_AGE = 3600        # seconds a plan stays valid for submission
MAX_SHOTS, MAX_CIRCUITS = 420000, 25
POLL = 10                  # seconds between two looks at a submitted job
REHEARSAL_DIVISOR = 10     # the rehearsal takes each circuit at a tenth of its shots


def pubs(rounds=ROUNDS, shots=SHOTS):
    # (label, round, shots) in the order the job runs them
    return [(label, r, shots[label]) for r in range(rounds) for label in H.LABELS]


def _seconds(stamp):
    # a UTC time stamp, with or without fractions of a second, as seconds since the epoch
    return calendar.timegm(time.strptime(stamp[:19], "%Y-%m-%dT%H:%M:%S"))


def safeties(plan, usage, now=None, record=None):
    # the reasons not to submit; none when the plan may run. record: the path a job record would take
    now = time.time() if now is None else now
    reasons = []
    if record is not None and os.path.exists(record):
        reasons.append(f"{os.path.basename(record)} exists; move the record of the earlier job before another run")
    age = now - _seconds(plan["made_utc"])
    if not -60 <= age <= PLAN_MAX_AGE:
        reasons.append(f"the plan is {age / 60:.0f} minutes old, the limit is {PLAN_MAX_AGE // 60}; run dry again")
    shots = sum(n for _label, _r, n in plan["pubs"])
    if ([list(x) for x in plan["pubs"]] != [list(x) for x in pubs()] or shots > MAX_SHOTS
            or len(plan["pubs"]) > MAX_CIRCUITS or len(plan["qubits"]) != QUBITS):
        reasons.append(f"the plan ({len(plan['pubs'])} circuits, {shots} shots, {len(plan['qubits'])} qubits) is not "
                       "this script's job")
    if plan["estimate_s"] > ESTIMATE_LIMIT:
        reasons.append(f"the estimate, {plan['estimate_s']:.0f} s, is above {ESTIMATE_LIMIT} s")
    if not plan["estimate_s"] <= plan["cap_s"] <= CAP_LIMIT:
        reasons.append(f"the cap, {plan['cap_s']} s, is not between the estimate and {CAP_LIMIT} s")
    left = usage.get("usage_remaining_seconds", 0)
    if usage.get("usage_limit_reached") or left - plan["cap_s"] < RESERVE:
        reasons.append(f"{left} s are left; a job that ran to its cap of {plan['cap_s']} s would leave less than "
                       f"{RESERVE} s")
    return reasons


def overrun(running_since, now, cap):
    # whether a job has been running for more than its cap plus 20 s
    return running_since is not None and now - running_since > cap + 20


def watch(job, cap):
    # follows a submitted job to its end; cancels it if the server's estimate or its running time passes the cap
    while True:
        status = str(job.status())
        if job.in_final_state():
            print(f"  job {job.job_id()}: {status}")
            return status
        try:
            estimate = (job.usage_estimation or {}).get("quantum_seconds")
        except Exception:
            estimate = None
        started = (job.metrics().get("timestamps") or {}).get("running")
        since = _seconds(started) if started else None
        if (estimate is not None and estimate > cap) or overrun(since, time.time(), cap):
            job.cancel()
            print(f"  job {job.job_id()} cancelled by the watch: server estimate {estimate} s, running since {started}, "
                  f"cap {cap} s")
            return "CANCELLED BY THE WATCH"
        print(f"  job {job.job_id()}: {status}, server estimate {estimate} s, running since {started}", flush=True)
        time.sleep(POLL)


def pooled(got, rounds):
    # {label: outcomes of rounds 0..rounds-1, one after another}
    return {label: np.concatenate([got[f"{label}_{r}"] for r in range(rounds)]) for label in H.LABELS}


def analyse(got, qubits, rounds, name, provenance):
    # the certificate from all rounds pooled, each qubit's coherence per round, and its bound after 1..rounds rounds
    from qrng_attest.quantum import attested_coherence_entropy, calibrate_readout, coherence_min_entropy
    total = pooled(got, rounds)
    cert, Ms = H.certify(total, name, provenance)
    att = cert.attestation
    nq = len(qubits)

    def bound(data, i):
        M, n = calibrate_readout(data["prep0"][:, i], data["prep1"][:, i])
        return attested_coherence_entropy(data["x"][:, i], data["y"][:, i], M, calib_shots=n)

    by_round = [[bound({label: got[f"{label}_{r}"] for label in H.LABELS}, i) for r in range(rounds)] for i in range(nq)]
    after = [[bound(pooled(got, r + 1), i) for r in range(rounds)] for i in range(nq)]
    lines = [f"  {'qubit':>5} {'F0':>7} {'F1':>7} {'P(0) Z':>7} {'coherence':>10} {'at least':>9} {'population':>11} "
             f"{'bound':>7} {'limit':>7}"]
    for q, M, row, z in zip(qubits, Ms, att.per_qubit, total["z"].T):
        lines.append(f"  {q:5d} {M[0, 0]:7.4f} {M[1, 1]:7.4f} {float(np.mean(z == 0)):7.4f} {row['coherence']:10.4f} "
                     f"{row['coherence_low']:9.4f} {row['h_population']:11.4f} {row['h_coherence']:7.4f} "
                     f"{coherence_min_entropy(row['coherence']):7.4f}")
    lines += ["", "  coherence measured in each round", f"  {'qubit':>5} " + " ".join(f"{r + 1:>7d}" for r in range(rounds))
              + f" {'spread':>7}"]
    for q, rows in zip(qubits, by_round):
        c = [x["coherence"] for x in rows]
        lines.append(f"  {q:5d} " + " ".join(f"{v:7.4f}" for v in c) + f" {max(c) - min(c):7.4f}")
    lines += ["", "  bound after the first rounds, pooled", f"  {'qubit':>5} " + " ".join(f"{r + 1:>7d}" for r in range(rounds))]
    for q, rows in zip(qubits, after):
        lines.append(f"  {q:5d} " + " ".join(f"{x['h_coherence']:7.4f}" for x in rows))
    lines += ["", str(cert)]
    data = {"attested": cert.rate, "passed": bool(cert.passed),
            "extractable": cert.plan.output_bits, "alpha": cert.plan.order,
            "per_qubit": [{k: row[k] for k in ("coherence", "coherence_low", "h_population", "h_coherence", "h_classical",
                                              "classical_track", "h_quantum", "s3_low", "s3_high", "h_long_run")}
                          for row in att.per_qubit],
            "readout": [[float(M[0, 0]), float(M[1, 1])] for M in Ms],
            "coherence_by_round": [[x["coherence"] for x in rows] for rows in by_round],
            "bound_after_rounds": [[x["h_coherence"] for x in rows] for rows in after]}
    return "\n".join(lines) + "\n", data


def dry(name=None, rehearse=True):
    from qiskit import qpy
    from qiskit.transpiler import generate_preset_pass_manager
    from qiskit_ibm_runtime import QiskitRuntimeService

    svc = QiskitRuntimeService()
    inst, usage = H.guards(svc, 1)
    print(f"instance {inst.get('name')}: plan {inst.get('plan')}, pricing {inst.get('pricing_type')}; "
          f"{usage['usage_remaining_seconds']} s of {usage['usage_limit_seconds']} left this period")
    backend = svc.backend(name) if name else svc.least_busy(operational=True, simulator=False)
    status, target = backend.status(), backend.target
    print(f"{backend.name}: operational {status.operational}, pending jobs {status.pending_jobs}, "
          f"{target.num_qubits} qubits")
    limit = backend.configuration().max_shots
    assert max(SHOTS.values()) <= limit, f"{max(SHOTS.values())} shots in one circuit, the backend takes {limit}"

    qubits = H.pick_qubits(backend, QUBITS)
    cal = H.calibration(target, qubits)
    print(f"qubits {qubits}")
    for q in qubits:
        c = cal[str(q)]
        print(f"  {q}: readout error {c['readout_error']:.4f}, sx error {c['sx_error']:.1e}, "
              f"T1 {c['t1'] * 1e6:.0f} us, T2 {c['t2'] * 1e6:.0f} us")

    pm = generate_preset_pass_manager(optimization_level=1, backend=backend, initial_layout=qubits,
                                      scheduling_method="alap")
    tqcs = [pm.run(c) for c in H.circuits(QUBITS)]
    for label, t in zip(H.LABELS, tqcs):
        assert H.measured(t) == qubits, f"circuit {label} reads {H.measured(t)}, not {qubits}"
    count = dict(zip(H.LABELS, (H.pulses(t) for t in tqcs)))
    assert count["x"] > count["z"] and count["y"] > count["z"], f"a basis change was cancelled: {count}"
    print(f"all {len(tqcs)} circuits on those qubits; sx and x pulses per circuit {count}")

    durs = dict(zip(H.LABELS, (t.estimate_duration(target, unit="s") for t in tqcs)))
    rep = backend.configuration().default_rep_delay
    rep = rep * 1e-6 if rep and rep > 1e-3 else (rep or 250e-6)
    plan_pubs = pubs()
    shots = sum(n for _l, _r, n in plan_pubs)
    est = sum(n * (durs[label] + rep) for label, _r, n in plan_pubs)
    need, cap = math.ceil(2 * est), math.ceil(1.5 * est) + 10
    assert cap <= CAP_LIMIT, f"the plan would cap the job at {cap} s, above {CAP_LIMIT} s"
    print(f"{ROUNDS} rounds of {dict(SHOTS)} shots: {len(plan_pubs)} circuits, {shots} shots; the backend takes "
          f"{limit} per circuit")
    print(f"scheduled durations {[round(durs[l] * 1e6, 2) for l in H.LABELS]} us; rep delay {rep * 1e6:.0f} us")
    print(f"estimate: {est:.0f} s of shots plus overhead; submit needs {need} s left and caps the job at {cap} s")
    H.guards(svc, need)

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "circuits.qpy"), "wb") as fh:
        qpy.dump(tqcs, fh)
    plan = {"backend": backend.name, "qubits": qubits, "rounds": ROUNDS, "shots_per_round": SHOTS,
            "pubs": [list(p) for p in plan_pubs], "labels": list(H.LABELS), "durations_s": durs, "rep_delay_s": rep,
            "estimate_s": est, "need_s": need, "cap_s": cap, "calibration": cal,
            "made_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(plan, open(os.path.join(OUT, "plan.json"), "w"), indent=1)

    try:
        if not rehearse:
            raise ImportError
        from qiskit_aer import AerSimulator
    except ImportError:
        print("no rehearsal" if not rehearse else "qiskit-aer is not installed; no rehearsal")
    else:
        print()
        print(f"rehearsal of the same circuits under the Aer noise model of {backend.name}, each at 1/{REHEARSAL_DIVISOR} "
              f"of its shots")
        sim = AerSimulator.from_backend(backend)
        got, seed = {}, 7
        for label, r, n in plan_pubs:
            n = n // REHEARSAL_DIVISOR
            mem = sim.run(tqcs[H.LABELS.index(label)], shots=n, memory=True, seed_simulator=seed).result().get_memory()
            seed += n + 1
            got[f"{label}_{r}"] = np.array([[int(x) for x in s.replace(" ", "")[::-1]] for s in mem], int)
        print(analyse(got, qubits, ROUNDS, f"{backend.name} model", "synthetic")[0])
    print(f"saved to {OUT}; nothing was submitted")


def submit(word=""):
    from qiskit import qpy
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    import qiskit
    import qiskit_ibm_runtime

    plan = json.load(open(os.path.join(OUT, "plan.json")))
    with open(os.path.join(OUT, "circuits.qpy"), "rb") as fh:
        tqcs = list(qpy.load(fh))
    need = math.ceil(2 * plan["estimate_s"])
    svc = QiskitRuntimeService()
    inst, usage = H.guards(svc, need)
    base = os.path.join(REPO, "results", f"hardware_coherence_rounds_{plan['backend']}")
    reasons = safeties(plan, usage, record=base + "_job.json")
    print(f"plan of {plan['made_utc']}: {plan['backend']}, qubits {plan['qubits']}, {len(plan['pubs'])} circuits, "
          f"{sum(n for _l, _r, n in plan['pubs'])} shots, estimate {plan['estimate_s']:.0f} s, cap {plan['cap_s']} s; "
          f"{usage['usage_remaining_seconds']} s left, plan {inst.get('plan')}, pricing {inst.get('pricing_type')}")
    for reason in reasons:
        print(f"  stop: {reason}")
    if reasons:
        raise SystemExit("not submitted")
    if word != "go":
        print("no reason to stop; not submitted without the word go")
        return
    backend = svc.backend(plan["backend"])
    if not backend.status().operational:
        raise SystemExit(f"{plan['backend']} is not operational; stopping")
    target = backend.target
    for q in plan["qubits"]:
        if target.qubit_properties[q] is None or (q,) not in target["measure"]:
            raise SystemExit(f"qubit {q} is no longer usable on {plan['backend']}; stopping")
    sampler = SamplerV2(mode=backend)
    sampler.options.max_execution_time = plan["cap_s"]
    sampler.options.dynamical_decoupling.enable = False
    sampler.options.twirling.enable_gates = False
    sampler.options.twirling.enable_measure = False
    cal = H.calibration(target, plan["qubits"])
    job_pubs = [(tqcs[plan["labels"].index(label)], None, n) for label, _r, n in plan["pubs"]]
    inst, usage = H.guards(svc, need)          # once more, immediately before the job
    if safeties(plan, usage, record=base + "_job.json"):
        raise SystemExit("not submitted")
    job = sampler.run(job_pubs)
    stub = {"job_id": job.job_id(), "backend": plan["backend"], "qubits": plan["qubits"], "rounds": plan["rounds"],
            "shots_per_round": plan["shots_per_round"], "pubs": plan["pubs"], "labels": plan["labels"],
            "submitted_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "plan": inst.get("plan"), "pricing_type": inst.get("pricing_type"), "estimate_s": plan["estimate_s"],
            "cap_s": plan["cap_s"], "calibration_at_submission": cal,
            "qiskit_version": qiskit.__version__, "qiskit_ibm_runtime_version": qiskit_ibm_runtime.__version__}
    json.dump(stub, open(base + "_job.json", "w"), indent=1)
    with open(base + "_circuits.qpy", "wb") as fh:
        qpy.dump(tqcs, fh)
    print(f"submitted job {stub['job_id']} to {plan['backend']}")
    stub["watch"] = watch(job, plan["cap_s"])
    json.dump(stub, open(base + "_job.json", "w"), indent=1)


def collect():
    from qiskit_ibm_runtime import QiskitRuntimeService

    plan = json.load(open(os.path.join(OUT, "plan.json")))
    base = os.path.join(REPO, "results", f"hardware_coherence_rounds_{plan['backend']}")
    stub = json.load(open(base + "_job.json"))
    if not os.path.exists(base + "_shots.npz"):
        job = QiskitRuntimeService().job(stub["job_id"])
        status = str(job.status())
        if "DONE" not in status.upper():
            print(f"job {stub['job_id']}: {status}; nothing collected")
            return
        got = {}
        for (label, r, n), pub in zip(stub["pubs"], job.result()):
            strings = next(iter(pub.data.values())).get_bitstrings()
            assert len(strings) == n, f"{label} of round {r}: {len(strings)} outcomes, {n} asked"
            got[f"{label}_{r}"] = np.array([[int(x) for x in s.replace(" ", "")[::-1]] for s in strings], np.uint8)
        np.savez_compressed(base + "_shots.npz", **got)
        stub["usage"] = job.metrics().get("usage", {})
        stub["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        json.dump(stub, open(base + "_job.json", "w"), indent=1)
    # everything below is computed from the saved shots, so it can be recomputed without the account
    got = {k: v.astype(int) for k, v in np.load(base + "_shots.npz").items()}
    body, data = analyse(got, stub["qubits"], stub["rounds"], stub["backend"], "measured")
    text = (f"job {stub['job_id']} on {stub['backend']}, qubits {stub['qubits']}, {stub['rounds']} rounds of "
            f"{stub['shots_per_round']} shots, plan {stub['plan']} and {stub['pricing_type']}, usage {stub.get('usage', {})}\n\n"
            + body)
    with open(base + ".txt", "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    json.dump({**stub, **data}, open(base + ".json", "w"), indent=1)
    print(text, end="")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "dry":
        names = [a for a in sys.argv[2:] if a != "norehearse"]
        dry(names[0] if names else None, rehearse="norehearse" not in sys.argv[2:])
    elif mode == "submit":
        submit(sys.argv[2] if len(sys.argv) > 2 else "")
    elif mode == "collect":
        collect()
    else:
        print("usage: python scripts/coherence_rounds.py dry [backend] [norehearse] | submit [go] | collect")
