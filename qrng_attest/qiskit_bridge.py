# sources, readout matrices and calibration from a Qiskit backend
import numpy as np

try:
    import qiskit
    _HAVE = True
except Exception:
    _HAVE = False


def available():
    return _HAVE


def _require():
    if not _HAVE:
        raise RuntimeError("qiskit_bridge needs qiskit (install in the project venv)")


def gate_duration(backend, name, qubits, default=None):
    # gate duration in seconds, from the BackendV2 Target first, then BackendV1 properties.
    try:
        inst = backend.target[name][tuple(qubits)]
        if inst is not None and inst.duration is not None:
            return float(inst.duration)
    except Exception:
        pass
    try:
        return float(backend.properties().gate_length(name, list(qubits)))
    except Exception:
        return default


def readout_matrix_from_backend(backend, qubit):
    # readout matrix from backend calibration. BackendV2 gives only a symmetric error; prefer self_calibrate_backend.
    _require()
    from .quantum import assignment_matrix
    try:
        p = backend.properties()
        F0 = 1.0 - float(p.qubit_property(qubit, "prob_meas1_prep0")[0])
        F1 = 1.0 - float(p.qubit_property(qubit, "prob_meas0_prep1")[0])
        return assignment_matrix(F0, F1)
    except Exception:
        pass
    try:
        err = float(backend.target["measure"][(qubit,)].error)
        return assignment_matrix(1.0 - err, 1.0 - err)
    except Exception:
        raise RuntimeError("could not read readout error for qubit " + str(qubit)
                           + "; use self_calibrate_backend so the certificate is not based on assumed perfect readout")


def self_calibrate_backend(backend, qubits=None, shots=4000):
    # measure the assignment matrix on the device itself from |0> and |1> preparations.
    _require()
    from qiskit import QuantumCircuit, transpile
    from .quantum import calibrate_readout
    qubits = list(range(backend.num_qubits)) if qubits is None else list(qubits)
    nq = len(qubits)
    outcomes = {}
    for prep in (0, 1):
        qc = QuantumCircuit(nq, nq)
        if prep == 1:
            for i in range(nq):
                qc.x(i)
        qc.measure(range(nq), range(nq))
        tqc = transpile(qc, backend, initial_layout=qubits)
        mem = _execute_shots(backend, tqc, shots)
        outcomes[prep] = np.array([[int(x) for x in m.replace(" ", "")[::-1]] for m in mem], int)
    Ms = []
    n_cal = shots
    for i in range(nq):
        M, n = calibrate_readout(outcomes[0][:, i], outcomes[1][:, i])
        Ms.append(M); n_cal = min(n_cal, n)
    return Ms, int(n_cal)


def source_circuit(nq, basis="z"):
    # |+> on every qubit, measured in Z, X or Y; the barrier keeps the transpiler from cancelling the preparation
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(nq, nq)
    for i in range(nq):
        qc.h(i)
    if basis != "z":
        qc.barrier()
        for i in range(nq):
            if basis == "y":
                qc.sdg(i)
            qc.h(i)
    qc.measure(range(nq), range(nq))
    return qc


def selftest_circuits(nq):
    # preparations 0, 1, + and - measured in Z and in X, keyed (x, y), for the dimension witness
    from qiskit import QuantumCircuit
    out = {}
    for x in range(4):
        for y in (0, 1):
            qc = QuantumCircuit(nq, nq)
            for q in range(nq):
                if x in (1, 3):
                    qc.x(q)
                if x in (2, 3):
                    qc.h(q)
            qc.barrier()
            if y == 1:
                qc.h(range(nq))
            qc.measure(range(nq), range(nq))
            out[(x, y)] = qc
    return out


def selftest_on_backend(backend, shots=4000, qubits=None, seed=0):
    # the eight circuits under the backend's noise model in Aer: per qubit, the witness and its entropy bound
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    from .quantum import self_testing_entropy
    qubits = list(range(2)) if qubits is None else list(qubits)
    sim = AerSimulator.from_backend(backend)
    ones = [[[0, 0] for _x in range(4)] for _q in qubits]
    for k, ((x, y), qc) in enumerate(selftest_circuits(len(qubits)).items()):
        tqc = transpile(qc, backend, initial_layout=qubits, optimization_level=1)
        memory = sim.run(tqc, shots=shots, memory=True, seed_simulator=seed + k * (shots + 1)).result().get_memory()
        for i in range(len(qubits)):
            ones[i][x][y] = sum(s.replace(" ", "")[::-1][i] == "1" for s in memory)
    return [self_testing_entropy(o, [[shots, shots]] * 4) for o in ones]


def _bits(memory):
    return np.array([[int(x) for x in s.replace(" ", "")[::-1]] for s in memory], int)


def quantum_source_from_backend(backend, n_shots=20000, qubits=None, seed=0, provenance="measured",
                                tomography=True):
    # run the |+> circuits under the backend's Aer noise model, returning a QuantumSource with readout matrices.
    _require()
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel
    from .quantum import QuantumArraySource

    qubits = list(range(backend.num_qubits)) if qubits is None else list(qubits)
    nq = len(qubits)
    sim = AerSimulator(noise_model=NoiseModel.from_backend(backend))
    got = {}
    for k, basis in enumerate(("z", "x", "y") if tomography else ("z",)):
        tqc = transpile(source_circuit(nq, basis), sim, initial_layout=qubits)
        # seeds further apart than the shot count: Aer seeds shot i of a run with seed + i
        got[basis] = _bits(sim.run(tqc, shots=n_shots, memory=True,
                                   seed_simulator=seed + k * (n_shots + 1)).result().get_memory())
    Ms = [readout_matrix_from_backend(backend, q) for q in qubits]
    name = getattr(backend, "name", "qiskit-backend")
    name = name() if callable(name) else name
    return QuantumArraySource(got["z"], readout=Ms, name=f"{name} QRNG", provenance=provenance,
                              x_shots=got.get("x"), y_shots=got.get("y"))


def _execute_shots(backend, tqc, n_shots):
    # get ordered per-shot bitstrings from the backend: try the runtime SamplerV2, then backend.run(memory=True)
    try:
        from qiskit_ibm_runtime import SamplerV2
        res = SamplerV2(mode=backend).run([tqc], shots=n_shots).result()
        bitarray = next(iter(res[0].data.values()))
        return bitarray.get_bitstrings()
    except Exception:
        return backend.run(tqc, shots=n_shots, memory=True).result().get_memory()


def run_on_backend(backend, n_shots=4096, qubits=None, provenance="measured", self_calibrate=True,
                   calibration_shots=4000, tomography=True):
    # run the |+> circuits ON the backend itself (real QPU or its runtime), returning a QuantumSource of
    _require()
    from qiskit import transpile
    from .quantum import QuantumArraySource
    qubits = list(range(backend.num_qubits)) if qubits is None else list(qubits)
    nq = len(qubits)
    n_cal = None
    if self_calibrate:
        Ms, n_cal = self_calibrate_backend(backend, qubits=qubits, shots=calibration_shots)
    else:
        Ms = [readout_matrix_from_backend(backend, q) for q in qubits]
    got = {}
    for basis in ("z", "x", "y") if tomography else ("z",):
        tqc = transpile(source_circuit(nq, basis), backend, initial_layout=qubits)
        got[basis] = _bits(_execute_shots(backend, tqc, n_shots))
    name = getattr(backend, "name", "qpu")
    name = name() if callable(name) else name
    if "fake" in type(backend).__module__.lower() or type(backend).__name__.lower().startswith("fake"):
        provenance = "synthetic"                          # a recorded device model, run locally
    return QuantumArraySource(got["z"], readout=Ms, name=f"{name} QRNG", provenance=provenance,
                              calibration_shots=n_cal, x_shots=got.get("x"), y_shots=got.get("y"))
