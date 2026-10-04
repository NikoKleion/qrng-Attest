# qrng_attest command line: certify a Qiskit backend, a capture file or the built-in demo source
import argparse
import os
import sys

import numpy as np


def _demo(args):
    from . import synthetic_qrng, synthetic_tomography, QuantumArraySource
    from .certify import certify_source
    shots, M = synthetic_qrng(args.shots, args.qubits, source_bias=args.bias,
                              crosstalk=args.crosstalk, drift=args.drift, seed=args.seed)
    x, y = (None, None) if args.no_tomography else synthetic_tomography(
        args.shots, args.qubits, source_bias=args.bias, coherence=args.coherence, seed=args.seed)
    src = QuantumArraySource(shots, readout=M, provenance="synthetic",
                             calibration_shots=args.calibration_shots, name="demo QRNG", x_shots=x, y_shots=y)
    print(certify_source(src))


def _backend(args):
    from .qiskit_bridge import available
    from .certify import certify_backend
    if not available():
        print("qiskit is required: pip install qiskit"); return 1
    backend = _load_backend(args.name)
    if backend is None:
        print(f"backend {args.name!r} is not among the fake backends of qiskit-ibm-runtime (for example manila, "
              f"fez, brisbane); a live device goes through certify_backend(backend) in Python")
        return 1
    qubits = None if args.qubits is None else [int(q) for q in args.qubits.split(",")]
    print(certify_backend(backend, shots=args.shots, qubits=qubits, simulate=args.simulate,
                          tomography=not args.no_tomography))


def _load_backend(name):
    try:
        import qiskit_ibm_runtime.fake_provider as fp
    except Exception:
        return None
    short = name[4:] if name.lower().startswith("ibm_") else name
    for cand in ("Fake" + short.capitalize() + "V2", "Fake" + short.capitalize(), "Fake" + short + "V2", "Fake" + short,
                 name):
        cls = getattr(fp, cand, None)
        if isinstance(cls, type):
            return cls()
    return None


def _capture(args):
    from . import assess, BytesFileSource
    if not os.path.isfile(args.path):
        print(f"no such file: {args.path}"); return 1
    print(assess(BytesFileSource(args.path, mode="bytes" if args.bytes else "bits"),
                 run_iid=not args.no_iid, iid_full=args.iid_full))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="qrng-attest",
                                 description="entropy assessment and device attestation for quantum RNG")
    sub = ap.add_subparsers(dest="cmd")

    d = sub.add_parser("demo", help="certify a built-in synthetic QRNG")
    d.add_argument("--shots", type=int, default=8000)
    d.add_argument("--qubits", type=int, default=4)
    d.add_argument("--bias", type=float, default=0.0)
    d.add_argument("--crosstalk", type=float, default=0.0)
    d.add_argument("--drift", type=float, default=0.0)
    d.add_argument("--calibration-shots", type=int, default=4000)
    d.add_argument("--coherence", type=float, default=1.0, help="coherence of the synthetic source state")
    d.add_argument("--no-tomography", action="store_true", help="leave out the X and Y basis shots")
    d.add_argument("--seed", type=int, default=0)
    d.set_defaults(fn=_demo)

    b = sub.add_parser("backend", help="certify a Qiskit backend by name")
    b.add_argument("name")
    b.add_argument("--shots", type=int, default=8192)
    b.add_argument("--simulate", action="store_true", help="run under the Aer noise model instead of the device")
    b.add_argument("--no-tomography", action="store_true", help="leave out the X and Y basis circuits")
    b.add_argument("--qubits", help="physical qubits, comma separated; the first five of the device when left out")
    b.set_defaults(fn=_backend)

    c = sub.add_parser("capture", help="assess a captured RNG file")
    c.add_argument("path")
    c.add_argument("--bytes", action="store_true", help="each byte is one sample; without it the file is a stream of bits")
    c.add_argument("--no-iid", action="store_true", help="leave out the IID track")
    c.add_argument("--iid-full", action="store_true", help="the 10000-shuffle permutation test of the IID track")
    c.set_defaults(fn=_capture)

    a = ap.parse_args(argv)
    if not getattr(a, "fn", None):
        ap.print_help(); return 0
    return a.fn(a) or 0


if __name__ == "__main__":
    sys.exit(main())
