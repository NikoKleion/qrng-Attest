# figures for the documentation, from the sweeps that write results/overcredit.txt and from a saved hardware record
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .overcredit import crosstalk_sweep, miscalibration_sweep, sweep

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "results", "figures")

INDIGO, TEAL, AMBER, SLATE = "#6366F1", "#14B8A6", "#F59E0B", "#94A3B8"
BLUE, NAVY = "#5170FF", "#004AAD"
ROUNDS = os.path.join(ROOT, "results", "hardware_coherence_rounds_ibm_fez.json")


def overcredit_figure(path, n_shots=60000):
    """The black-box credit against the attested bound, and the cost of mis-calibration."""
    points = sweep(n_shots=n_shots)
    miscal = miscalibration_sweep(n_shots=n_shots)

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.1))

    F = [p["F"] for p in points]
    truth = points[0]["h_true"]
    credited = [p["credited"] for p in points]
    attested = [p["attested"] for p in points]
    ax.fill_between(F, [truth] * len(F), credited, color=AMBER, alpha=0.18, lw=0, label="over-credit")
    ax.plot(F, credited, "o-", color=AMBER, lw=2.2, ms=6, label="SP 800-90B, output only")
    ax.plot(F, attested, "s-", color=INDIGO, lw=2.2, ms=6, label="attested through the calibration")
    ax.axhline(truth, color=SLATE, ls="--", lw=1.8, label=f"true source entropy, {truth:.4f}")
    ax.set_xlabel("readout fidelity")
    ax.set_ylabel("min-entropy (bits per bit)")
    ax.set_title("a near-deterministic source read through a noisy detector", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6)

    dF = [m["dF"] for m in miscal]
    bx.plot(dF, [m["attested"] for m in miscal], "o-", color=AMBER, lw=2.2, ms=6,
            label="calibration trusted as exact")
    bx.plot(dF, [m["attested_boxed"] for m in miscal], "s-", color=INDIGO, lw=2.2, ms=6,
            label="calibration margin propagated")
    bx.axhline(miscal[0]["h_true"], color=SLATE, ls="--", lw=1.8, label="true source entropy")
    bx.set_xlabel("error in the assumed readout fidelity")
    bx.set_ylabel("attested min-entropy (bits per bit)")
    bx.set_title("sensitivity to an overstated calibration", fontsize=11)
    bx.set_ylim(bottom=0)
    bx.legend(fontsize=8.5, frameon=False)
    bx.grid(alpha=0.25, lw=0.6)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return points, miscal


def crosstalk_figure(path, n_shots=60000):
    """Per-stream assessment against the joint rate when two qubits are correlated."""
    points = crosstalk_sweep(n_shots=n_shots)
    c = [p["c"] for p in points]

    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot(c, [p["credited_rate"] for p in points], "o-", color=AMBER, lw=2.2, ms=6,
            label="per stream, assessed separately")
    ax.plot(c, [p["mi_debit_rate"] for p in points], "^-", color=SLATE, lw=2, ms=6,
            label="one minus the mutual information")
    ax.plot(c, [p["joint_bound_rate"] for p in points], "s-", color=INDIGO, lw=2.2, ms=6,
            label="worst pair joint bound")
    ax.plot(c, [p["truth_rate"] for p in points], "--", color="#DB2777", lw=2, label="true joint rate")
    ax.set_xlabel("probability that the second qubit copies the first")
    ax.set_ylabel("rate (bits per shot)")
    ax.set_title("two correlated qubits", fontsize=11)
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return points


def hardware_figure(path, record=ROUNDS):
    """The rounds job: each qubit's bound as its rounds are pooled, and what the pooled record certifies."""
    meta = json.load(open(record))
    per_round = meta["shots_per_round"]["x"]
    shots = [per_round * (r + 1) for r in range(meta["rounds"])]
    after = meta["bound_after_rounds"]
    qubits = meta["qubits"]
    one = [q["h_coherence"] for q in meta["per_qubit"]]
    long_run = [q["h_long_run"] for q in meta["per_qubit"]]
    order = sorted(range(len(qubits)), key=lambda i: one[i])

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.1), gridspec_kw={"width_ratios": [1, 1.25]})
    for i in range(len(qubits)):
        ax.plot(shots, after[i], "-", color=BLUE, lw=1.1, alpha=0.55)
    ax.plot(shots, after[order[0]], "o-", color=NAVY, lw=2.2, ms=5, label=f"qubit {qubits[order[0]]}, the least")
    ax.plot(shots, after[order[-1]], "s-", color="black", lw=2.2, ms=5, label=f"qubit {qubits[order[-1]]}, the largest")
    ax.set_xlabel("shots per circuit, rounds pooled")
    ax.set_ylabel("coherence bound, one shot (bits per bit)")
    ax.set_title(f"{meta['backend']}, {len(qubits)} qubits, {meta['rounds']} rounds", fontsize=11)
    ax.set_xticks(shots)
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6)

    x = range(len(qubits))
    bx.bar(x, [long_run[i] for i in order], color=BLUE, alpha=0.35, width=0.72, label="long outputs")
    bx.bar(x, [one[i] for i in order], color=NAVY, width=0.72, label="one shot")
    rate = meta["extractable"] / (len(qubits) * meta["shots_per_round"]["z"] * meta["rounds"])
    bx.axhline(rate, color="black", ls="--", lw=1.6, label=f"over this output, {rate:.4f}")
    bx.set_xticks(list(x))
    bx.set_xticklabels([str(qubits[i]) for i in order], fontsize=8.5)
    bx.set_xlabel("qubit")
    bx.set_ylabel("certified (bits per bit)")
    bx.set_ylim(0.7, 1.04)
    bx.set_title("all rounds pooled", fontsize=11)
    bx.legend(fontsize=8.5, frameon=False, loc="upper center", ncol=3)
    bx.grid(alpha=0.25, lw=0.6, axis="y")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return {"shots": shots, "after": after, "one_shot": one, "long_run": long_run, "rate": rate, "qubits": qubits}


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    drawn = hardware_figure(os.path.join(FIGDIR, "fig3_hardware_rounds.png"))
    print("fig3_hardware_rounds.png")
    print(f"  one shot {min(drawn['one_shot']):.4f} to {max(drawn['one_shot']):.4f}, long outputs "
          f"{min(drawn['long_run']):.4f} to {max(drawn['long_run']):.4f}, over the output {drawn['rate']:.4f}")
    points, miscal = overcredit_figure(os.path.join(FIGDIR, "fig1_overcredit.png"))
    print("fig1_overcredit.png")
    for p in points:
        print(f"  F={p['F']:.3f}  credited {p['credited']:.4f}  attested {p['attested']:.4f}")
    for m in miscal:
        print(f"  dF={m['dF']:+.2f}  point {m['attested']:.4f}  boxed {m['attested_boxed']:.4f}")
    cross = crosstalk_figure(os.path.join(FIGDIR, "fig2_crosstalk.png"))
    print("fig2_crosstalk.png")
    for p in cross:
        print(f"  c={p['c']:.2f}  per stream {p['credited_rate']:.4f}  "
              f"joint {p['joint_bound_rate']:.4f}  truth {p['truth_rate']:.4f}")


if __name__ == "__main__":
    main()
