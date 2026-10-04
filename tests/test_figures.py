# the figures come from the same sweeps as results/overcredit.txt and agree with it
import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# results/overcredit.txt, readout fidelity against the SP 800-90B estimate and the attested bound
OVERCREDIT = {0.55: (0.8580, 0.0000), 0.75: (0.4247, 0.0145),
              0.95: (0.0974, 0.0247), 0.999: (0.0282, 0.0268)}
TRUE_ENTROPY = 0.0291


def _figures():
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        raise unittest.SkipTest("requires matplotlib")
    from qrng_attest import figures
    return figures


class FigureData(unittest.TestCase):
    def test_overcredit_sweep_matches_the_saved_table(self):
        _figures()
        from qrng_attest.overcredit import sweep
        for point in sweep(n_shots=60000):
            want = OVERCREDIT.get(round(point["F"], 3))
            if want is None:
                continue
            self.assertAlmostEqual(point["credited"], want[0], places=3)
            self.assertAlmostEqual(point["attested"], want[1], places=3)
            self.assertLessEqual(point["attested"], TRUE_ENTROPY + 1e-9,
                                 f"attested above the truth at F={point['F']}")

    def test_the_joint_bound_never_exceeds_the_true_rate(self):
        _figures()
        from qrng_attest.overcredit import crosstalk_sweep
        for point in crosstalk_sweep(n_shots=60000):
            self.assertLessEqual(point["joint_bound_rate"], point["truth_rate"] + 1e-9,
                                 f"joint bound above the truth at c={point['c']}")

    def test_hardware_figure_draws_the_saved_record(self):
        figures = _figures()
        import json
        import tempfile
        meta = json.load(open(figures.ROUNDS))
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "hardware.png")
            drawn = figures.hardware_figure(path)
            self.assertGreater(os.path.getsize(path), 10_000)
        self.assertEqual(drawn["shots"], [20000, 40000, 60000, 80000, 100000])
        self.assertEqual(drawn["qubits"], meta["qubits"])
        self.assertEqual(drawn["after"], meta["bound_after_rounds"])
        for i, row in enumerate(meta["per_qubit"]):
            self.assertEqual(drawn["one_shot"][i], row["h_coherence"])
            self.assertAlmostEqual(drawn["after"][i][-1], row["h_coherence"], places=12)
            self.assertLess(row["h_coherence"], row["h_long_run"])
        self.assertAlmostEqual(drawn["rate"], meta["extractable"] / 240000, places=12)

    def test_the_readme_figures_are_present(self):
        figdir = os.path.join(ROOT, "results", "figures")
        for name in ("fig1_overcredit.png", "fig2_crosstalk.png", "fig3_hardware_rounds.png"):
            path = os.path.join(figdir, name)
            self.assertTrue(os.path.exists(path), path)
            self.assertGreater(os.path.getsize(path), 10_000, path)


if __name__ == "__main__":
    unittest.main()
