# submit the QRNG job to an IBM QPU and return without waiting
import sys

from qiskit import transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

from run_hardware import build_circuits, QUBITS, SHOTS


def main():
    svc = QiskitRuntimeService()
    inst = svc.instances()[0]
    if inst.get("pricing_type") != "free" or inst.get("plan") != "open":
        print("instance is not the free open plan, stopping:", inst.get("plan"), inst.get("pricing_type"))
        return 1
    be = svc.least_busy(operational=True, simulator=False)
    circs = transpile(build_circuits(len(QUBITS)), be, initial_layout=QUBITS)
    job = SamplerV2(mode=be).run(circs, shots=SHOTS)
    jid = job.job_id()
    with open("last_job.txt", "w", encoding="utf-8") as f:
        f.write(be.name + "\n" + jid + "\n")
    print("submitted to", be.name)
    print("job id", jid)
    print("run scripts/poll_job.py to check status and get the certificate.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
