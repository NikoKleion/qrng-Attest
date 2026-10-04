# the self-testing bound: the dimension witness of eight prepare-and-measure settings and its min-entropy
import math

import numpy as np

from _harness import skip
from qrng_attest import quantum as Q


def _unit(rng):
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)


def _strategy(rng, parts):
    # independent devices: a mixture of preparation devices (four vectors each) and one of measurement devices
    q, r = rng.dirichlet(np.ones(parts)), rng.dirichlet(np.ones(parts))
    return (q, [[rng.random() * _unit(rng) for _ in range(4)] for _ in range(parts)],
            r, [[_unit(rng) for _ in range(2)] for _ in range(parts)])


def _witness_and_guess(q, S, r, T):
    # the witness of the observed probabilities, and the guessing probability of someone who knows both internal states
    p1 = np.zeros((4, 2))
    guess = 0.0
    for qk, Sk in zip(q, S):
        for rk, Tk in zip(r, T):
            for x in range(4):
                for y in range(2):
                    e = float(np.dot(Sk[x], Tk[y]))
                    p1[x, y] += qk * rk * (1.0 - e) / 2.0
                    guess += qk * rk * (1.0 + abs(e)) / 2.0 / 8.0
    W = abs((p1[0, 0] - p1[1, 0]) * (p1[2, 1] - p1[3, 1]) - (p1[2, 0] - p1[3, 0]) * (p1[0, 1] - p1[1, 1]))
    return W, guess


def test_witness_bound_values():
    # Lunghi et al. 2015: W = 1 certifies a min-entropy of about 0.2284, and W = 0 certifies nothing
    assert abs(-math.log2(Q.witness_guess(1.0)) - 0.2284) < 5e-5 and Q.witness_guess(0.0) == 1.0
    values = [Q.witness_guess(w) for w in np.linspace(0.0, 1.0, 50)]
    assert all(b <= a for a, b in zip(values, values[1:]))


def test_no_qubit_strategy_guesses_better_than_the_bound():
    # mixtures of qubit preparations and of measurements, drawn independently: guessing is at or below eq. 5
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(3000):
        W, guess = _witness_and_guess(*_strategy(rng, int(rng.integers(1, 4))))
        assert 0.0 <= W <= 1.0 + 1e-12
        worst = max(worst, guess - Q.witness_guess(W))
    assert worst <= 1e-12, worst
    # the BB84 settings reach W = 1 and the bound exactly
    z, x = np.array([0, 0, 1.0]), np.array([1.0, 0, 0])
    W, guess = _witness_and_guess([1.0], [[z, -z, x, -x]], [1.0], [[z, x]])
    assert abs(W - 1.0) < 1e-12 and guess <= Q.witness_guess(1.0) + 1e-12 and abs(guess - 0.75) < 1e-12


def test_classical_preparations_certify_nothing():
    # preparations diagonal in one basis give W = 0, whatever the outcomes look like
    z = np.array([0, 0, 1.0])
    rng = np.random.default_rng(2)
    for _ in range(20):
        S = [[(2 * rng.random() - 1) * z for _ in range(4)]]
        W, _guess = _witness_and_guess([1.0], S, [1.0], [[_unit(rng), _unit(rng)]])
        assert W < 1e-12
    n = 50000
    coin = [[n // 2, n // 2]] * 4
    r = Q.self_testing_entropy(coin, [[n, n]] * 4)
    assert r["witness_low"] == 0.0 and r["h_min"] == 0.0


def test_lower_limit_on_the_witness_holds_and_moves_with_the_counts():
    rng = np.random.default_rng(3)
    p = np.array([[0.01, 0.5], [0.98, 0.5], [0.5, 0.015], [0.5, 0.985]])
    true = abs((p[0, 0] - p[1, 0]) * (p[2, 1] - p[3, 1]) - (p[2, 0] - p[3, 0]) * (p[0, 1] - p[1, 1]))
    above = 0
    for _ in range(300):
        ones = rng.binomial(4000, p)
        r = Q.self_testing_entropy(ones.tolist(), [[4000, 4000]] * 4)
        above += r["witness_low"] > true
        assert r["witness_low"] <= r["witness"] <= 1.0
    assert above <= 3
    small = Q.self_testing_entropy(rng.binomial(1000, p).tolist(), [[1000, 1000]] * 4)["h_min"]
    large = Q.self_testing_entropy(rng.binomial(100000, p).tolist(), [[100000, 100000]] * 4)["h_min"]
    assert small < large < -math.log2(Q.witness_guess(true)) + 1e-9


def test_circuits_run_under_a_backend_noise_model():
    try:
        import qiskit_aer  # noqa: F401
        from qiskit_ibm_runtime.fake_provider import FakeManilaV2
    except ImportError:
        skip("requires qiskit-aer and qiskit-ibm-runtime")
    from qrng_attest.qiskit_bridge import selftest_circuits, selftest_on_backend
    circuits = selftest_circuits(2)
    assert sorted(circuits) == [(x, y) for x in range(4) for y in (0, 1)]
    assert circuits[(3, 1)].count_ops().get("h", 0) == 4 and circuits[(0, 0)].count_ops().get("h", 0) == 0
    out = selftest_on_backend(FakeManilaV2(), shots=6000, qubits=[0, 1], seed=1)
    assert len(out) == 2 and all(0.7 < r["witness_low"] <= r["witness"] < 1.0 and 0.04 < r["h_min"] < 0.2285 for r in out)
    assert all(r["p"][0][0] < 0.1 and r["p"][1][0] > 0.9 and r["p"][2][1] < 0.1 and r["p"][3][1] > 0.9 for r in out)
