# quantum-device min-entropy from Bloch geometry and readout deconvolution
import math
import numpy as np

from .sources import Source


# Bloch-sphere geometry
def bloch_z(p0):
    # Bloch z-component from P(0): r_z = 2 P0 - 1
    return 2.0 * float(p0) - 1.0


def geometric_min_entropy(p0):
    # min-entropy of a Z-basis outcome: -log2(max(P0, P1)).
    p0 = min(max(float(p0), 0.0), 1.0)
    return float(-math.log2(max(p0, 1.0 - p0, 1e-12)) + 0.0)


def state_min_entropy(theta):
    # min-entropy for a state at polar angle theta: P0 = cos^2(theta/2).
    return geometric_min_entropy(math.cos(theta / 2.0) ** 2)


# readout channel
def assignment_matrix(F0, F1):
    # assignment matrix M[a, b] = P(read a | true b).
    return np.array([[F0, 1.0 - F1], [1.0 - F0, F1]], float)


def readout_deconvolve(p_obs, M):
    # recover the true distribution from the observed: p_true = M^-1 p_obs, projected onto the simplex.
    p_obs = np.asarray(p_obs, float)
    det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    if abs(det) < 1e-9:
        return p_obs
    p_true = np.linalg.solve(M, p_obs)
    p_true = np.clip(p_true, 0.0, 1.0)
    s = p_true.sum()
    return p_true / s if s > 0 else np.array([0.5, 0.5])


def calibrate_readout(prep0_outcomes, prep1_outcomes):
    # measure the assignment matrix from known preparations instead of trusting a vendor number.
    a0 = np.asarray(prep0_outcomes, int).ravel()
    a1 = np.asarray(prep1_outcomes, int).ravel()
    n0 = max(len(a0), 1); n1 = max(len(a1), 1)
    F0 = float(np.count_nonzero(a0 == 0)) / n0
    F1 = float(np.count_nonzero(a1 == 1)) / n1
    return assignment_matrix(F0, F1), int(min(n0, n1))


Z_CONF = 2.5758293035489004


def _conservative_p0(p0, sigma, z=Z_CONF):
    # push bias away from 1/2 by z sigma for a one-sided lower bound.
    if p0 >= 0.5:
        return min(1.0, p0 + z * sigma)
    return max(0.0, p0 - z * sigma)


def _corner_entropy(p0_obs, sigma_obs, F0, F1):
    # deconvolved conservative entropy for one candidate assignment matrix.
    M = assignment_matrix(F0, F1)
    p_true = readout_deconvolve(np.array([p0_obs, 1.0 - p0_obs]), M)
    det = F0 * F1 - (1.0 - F0) * (1.0 - F1)
    sigma_src = sigma_obs / abs(det) if abs(det) > 1e-9 else sigma_obs
    p0_src = float(p_true[0])
    return geometric_min_entropy(_conservative_p0(p0_src, sigma_src)), p0_src


def attested_qubit_entropy(bits, M=None, calib_shots=None, M_margin=None):
    # quantum min-entropy of one qubit's stream, 99% one-sided lower bound after readout deconvolution.
    bits = np.asarray(bits, int)
    n = max(len(bits), 2)
    p0_obs = float(np.mean(bits == 0))
    h_obs = geometric_min_entropy(p0_obs)
    sigma_obs = math.sqrt(max(p0_obs * (1.0 - p0_obs), 1e-12) / (n - 1))
    if M is None:
        p0_cons = _conservative_p0(p0_obs, sigma_obs)
        h_q = geometric_min_entropy(p0_cons)
        return {"p0_obs": p0_obs, "h_observed": h_obs, "p0_source": p0_obs, "h_quantum": h_q,
                "r_z": bloch_z(p0_obs), "readout_dominated": False, "has_readout_model": False}
    F0, F1 = float(M[0, 0]), float(M[1, 1])
    if M_margin is not None:
        d0 = d1 = float(M_margin)
    elif calib_shots is not None:
        nc = max(int(calib_shots), 2)
        d0 = Z_CONF * math.sqrt(max(F0 * (1.0 - F0), 1e-12) / nc)
        d1 = Z_CONF * math.sqrt(max(F1 * (1.0 - F1), 1e-12) / nc)
    else:
        d0 = d1 = 0.0
    h_q, p0_src = _corner_entropy(p0_obs, sigma_obs, F0, F1)
    if d0 > 0.0 or d1 > 0.0:
        for f0 in (F0 - d0, F0, F0 + d0):
            for f1 in (F1 - d1, F1, F1 + d1):
                f0c = min(0.9995, max(0.505, f0)); f1c = min(0.9995, max(0.505, f1))
                h_c, _p = _corner_entropy(p0_obs, sigma_obs, f0c, f1c)
                if h_c < h_q:
                    h_q = h_c
    dominated = (h_obs - h_q) > 0.1
    return {"p0_obs": p0_obs, "h_observed": h_obs, "p0_source": p0_src, "h_quantum": h_q,
            "r_z": bloch_z(p0_src), "readout_dominated": dominated, "has_readout_model": True,
            "calib_margin": (d0, d1)}


def coherence_min_entropy(c):
    # min-entropy of a Z outcome for a state of coherence c = |S1 - i S2| (Fiorentino et al. 2007, eq. 7)
    c = min(max(float(c), 0.0), 1.0)
    return float(-math.log2((1.0 + math.sqrt(1.0 - c * c)) / 2.0) + 0.0)


def _calibration_margins(F0, F1, calib_shots=None, M_margin=None):
    # half-widths of the confidence box on the two readout fidelities
    if M_margin is not None:
        return float(M_margin), float(M_margin)
    if calib_shots is None:
        return 0.0, 0.0
    nc = max(int(calib_shots), 2)
    return (Z_CONF * math.sqrt(max(F0 * (1.0 - F0), 1e-12) / nc),
            Z_CONF * math.sqrt(max(F1 * (1.0 - F1), 1e-12) / nc))


def _wilson(p, n, z=Z_CONF):
    # Wilson score interval for a proportion; it stays inside (0, 1) when p is 0 or 1
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def _stokes_corners(bits, M=None, calib_shots=None, M_margin=None):
    # readout-corrected S = 2 P(0) - 1 of one basis, and its range at each corner of the calibration box
    bits = np.asarray(bits, int)
    n = max(len(bits), 2)
    p = float(np.mean(bits == 0))
    lo, hi = _wilson(p, n)
    F0, F1 = (1.0, 1.0) if M is None else (float(M[0, 0]), float(M[1, 1]))
    d0, d1 = _calibration_margins(F0, F1, calib_shots, M_margin)
    component = lambda q, f0, f1: (2.0 * q - 1.0 - (f0 - f1)) / (f0 + f1 - 1.0)
    # up to 1: an error seen in the calibration may be a fault of the prepared state, not of the readout
    corners = [(component(lo, f0, f1), component(hi, f0, f1))
               for f0 in {min(1.0, max(0.505, F0 - d0)), 1.0}
               for f1 in {min(1.0, max(0.505, F1 - d1)), 1.0}]
    centre = component(p, min(1.0, max(0.505, F0)), min(1.0, max(0.505, F1)))
    return float(min(1.0, max(-1.0, centre))), corners


def stokes_lower(bits, M=None, calib_shots=None, M_margin=None):
    # readout-corrected Stokes component S = 2 P(0) - 1 of one basis, and a lower bound on |S|
    centre, corners = _stokes_corners(bits, M, calib_shots, M_margin)
    if any(a <= 0.0 <= b for a, b in corners) or len({a > 0.0 for a, _b in corners}) > 1:
        return centre, 0.0
    return centre, float(min(1.0, min(min(abs(a), abs(b)) for a, b in corners)))


def stokes_upper(bits, M=None, calib_shots=None, M_margin=None):
    # an upper bound on |S| of one basis, over the same interval and corners
    _centre, corners = _stokes_corners(bits, M, calib_shots, M_margin)
    return float(min(1.0, max(max(abs(a), abs(b)) for a, b in corners)))


def attested_coherence_entropy(x_bits, y_bits=None, M=None, calib_shots=None, M_margin=None):
    # lower bound on the min-entropy of the Z outcome from X and Y basis shots of the same state
    s1, l1 = stokes_lower(x_bits, M, calib_shots, M_margin)
    s2, l2 = (0.0, 0.0) if y_bits is None else stokes_lower(y_bits, M, calib_shots, M_margin)
    c, c_low = min(1.0, math.hypot(s1, s2)), min(1.0, math.hypot(l1, l2))
    return {"s1": s1, "s2": s2, "coherence": c, "coherence_low": c_low,
            "h_coherence": coherence_min_entropy(c_low), "bases": "x" if y_bits is None else "xy"}


def _h2(x):
    return 0.0 if x <= 0.0 or x >= 1.0 else -x * math.log2(x) - (1.0 - x) * math.log2(1.0 - x)


def renyi_entropy(c, s3, alpha):
    # alpha entropy of a Z outcome given a purification, at coherence c and S3 (Tomamichel et al. 2009, eq. 16)
    s3 = min(abs(float(s3)), 1.0)
    c = min(max(float(c), 0.0), math.sqrt(1.0 - s3 * s3))
    p0, p1 = (1.0 + s3) / 2.0, (1.0 - s3) / 2.0
    if p1 <= 0.0:
        return 0.0
    limit = _h2(p0) - _h2((1.0 + math.sqrt(min(1.0, c * c + s3 * s3))) / 2.0)
    if alpha == 1.0:
        return float(limit + 0.0)
    t = p0 ** alpha + p1 ** alpha
    d = (p0 * p1) ** alpha * max(0.0, 1.0 - c * c / (1.0 - s3 * s3))
    top = (t + math.sqrt(max(0.0, t * t - 4.0 * d))) / 2.0
    value = alpha / (1.0 - alpha) * math.log2(top ** (1.0 / alpha) + (d / top) ** (1.0 / alpha))
    return float(min(value, limit) + 0.0)      # the entropy falls with alpha; the limit guards the rounding near 1


ORDERS = 1.0 + np.logspace(-5.0, 0.0, 501)      # orders alpha in (1, 2] tried for the bound on n shots


def worst_entropy(c_low, s_low, s_high, alpha):
    # least alpha entropy over c >= c_low and s_low <= |S3| <= s_high; in |S3| it rises and then falls
    return min(renyi_entropy(c_low, s_low, alpha), renyi_entropy(c_low, s_high, alpha))


def shot_entropy(qubits, n, epsilon):
    # lower bound on the epsilon-smooth min-entropy of n shots of each qubit (c_low, s_low, s_high): (bits, alpha)
    cost = 1.0 - 2.0 * math.log2(epsilon)
    best, order = 0.0, None
    for a in ORDERS:
        k = n * sum(worst_entropy(c, lo, hi, a) for c, lo, hi in qubits) - cost / (a - 1.0)
        if k > best:
            best, order = k, float(a)
    return float(best), order


def witness_guess(w):
    # largest average guessing probability at dimension witness W (Lunghi et al. 2015, eq. 5)
    w = min(max(float(w), 0.0), 1.0)
    return 0.5 * (1.0 + math.sqrt((1.0 + math.sqrt(1.0 - w * w)) / 2.0))


def self_testing_entropy(ones, shots):
    # ones[x][y] outcomes 1 in shots[x][y] rounds of preparation x in 0..3 and measurement y in 0..1 (eq. 4)
    p = [[ones[x][y] / shots[x][y] for y in (0, 1)] for x in range(4)]
    box = [[_wilson(p[x][y], shots[x][y]) for y in (0, 1)] for x in range(4)]
    det = lambda a, b, c, d: a * d - b * c
    centre = abs(det(p[0][0] - p[1][0], p[2][0] - p[3][0], p[0][1] - p[1][1], p[2][1] - p[3][1]))
    # each entry of the determinant is a difference of two proportions; the determinant is extreme at corners
    span = lambda x, x2, y: (box[x][y][0] - box[x2][y][1], box[x][y][1] - box[x2][y][0])
    A, B, C, D = span(0, 1, 0), span(2, 3, 0), span(0, 1, 1), span(2, 3, 1)
    corners = [det(a, b, c, d) for a in A for b in B for c in C for d in D]
    low = 0.0 if min(corners) <= 0.0 <= max(corners) else min(abs(min(corners)), abs(max(corners)))
    return {"witness": float(min(1.0, centre)), "witness_low": float(min(1.0, low)),
            "h_min": float(-math.log2(witness_guess(low)) + 0.0), "p": p}


# QuantumSource
class QuantumSource:
    name = "quantum-source"
    provenance = "unknown"
    expected_state = "plus"
    calibration_shots = None

    def shots(self):
        # (n_shots, n_qubits) array of measurement outcomes in {0, 1}
        raise NotImplementedError

    @property
    def n_qubits(self):
        return np.asarray(self.shots()).shape[1]

    def readout_matrices(self):
        # per-qubit 2x2 assignment matrices, or None.
        return getattr(self, "_readout", None)

    def basis_shots(self):
        # {"x": shots, "y": shots or None} of the same state measured in the X and Y bases, or None
        return getattr(self, "_basis", None)

    def qubit_stream(self, q):
        # one qubit's bit stream.
        return np.asarray(self.shots())[:, q]

    def describe(self):
        s = np.asarray(self.shots())
        flags = []
        if self.provenance == "synthetic":
            flags.append("synthetic quantum source: this run is simulated, not measured hardware")
        if self.basis_shots() is None:
            flags.append("no X or Y basis shots: the value assumes a pure source state and is an upper bound "
                         "otherwise")
        if self.readout_matrices() is None:
            flags.append("no readout model provided; quantum min-entropy falls back to the observed bias")
        elif self.calibration_shots is None:
            flags.append("readout matrix trusted as exact; supply calibration_shots to bound calibration error")
        return {"name": self.name, "n_shots": int(s.shape[0]), "n_qubits": int(s.shape[1]),
                "provenance": self.provenance, "expected_state": self.expected_state, "flags": flags}


class QuantumArraySource(QuantumSource):
    def __init__(self, shots, readout=None, name="qubit-readout", provenance="measured", expected_state="plus",
                 calibration_shots=None, x_shots=None, y_shots=None):
        self._s = np.ascontiguousarray(np.asarray(shots)).astype(int)
        if self._s.ndim == 1:
            self._s = self._s[:, None]
        nq = self._s.shape[1]
        if readout is not None and np.shape(readout) == (2, 2):
            readout = [np.asarray(readout, float)] * nq          # one matrix for every qubit
        if readout is not None and len(readout) != nq:
            raise ValueError(f"readout holds {len(readout)} matrices for {nq} qubits; give one 2x2 matrix or one per qubit")
        self._readout = readout
        self._basis = None
        if x_shots is not None:
            two = lambda a: None if a is None else np.asarray(a, int).reshape(len(a), -1)
            self._basis = {"x": two(x_shots), "y": two(y_shots)}
            for key, value in self._basis.items():
                if value is not None and value.shape[1] != nq:
                    raise ValueError(f"{key}_shots has {value.shape[1]} columns for {nq} qubits; shape is (shots, qubits)")
        self.name = name; self.provenance = provenance; self.expected_state = expected_state
        self.calibration_shots = calibration_shots

    def shots(self):
        return self._s


class QuantumCallableSource(QuantumSource):
    # live device: fn(n_shots) -> (n_shots, n_qubits) outcomes.
    def __init__(self, fn, n_shots, readout=None, name="live-qpu", provenance="measured", expected_state="plus"):
        self.fn = fn; self.n = int(n_shots); self._readout = readout
        self.name = name; self.provenance = provenance; self.expected_state = expected_state
        self._cache = None

    def shots(self):
        if self._cache is None:
            s = np.asarray(self.fn(self.n)).astype(int)
            self._cache = s[:, None] if s.ndim == 1 else s
        return self._cache


def synthetic_qrng(n_shots, n_qubits=4, source_bias=0.0, readout_F=(0.99, 0.98), crosstalk=0.0,
                   drift=0.0, seed=0, model_readout=True):
    # QRNG device model for testing and demos; returns (shots, readout_matrices or None).
    rng = np.random.default_rng(seed)
    F0, F1 = readout_F
    M = assignment_matrix(F0, F1)
    shots = np.empty((n_shots, n_qubits), int)
    ramp = np.linspace(0.0, drift, n_shots)
    for t in range(n_shots):
        b = source_bias + ramp[t]
        p0 = np.clip(0.5 + b, 0.02, 0.98)
        true = (rng.random(n_qubits) > p0).astype(int)
        if crosstalk > 0 and n_qubits > 1:
            copy = rng.random(n_qubits) < crosstalk
            true[1:] = np.where(copy[1:], true[:-1], true[1:])
        obs = np.where(rng.random(n_qubits) < np.where(true == 0, 1 - F0, 1 - F1), 1 - true, true)
        shots[t] = obs
    return shots, ([M] * n_qubits if model_readout else None)


def synthetic_tomography(n_shots, n_qubits=4, source_bias=0.0, coherence=1.0, readout_F=(0.99, 0.98), seed=0):
    # X and Y basis shots of the synthetic_qrng state: S3 = 2 bias, S1 = coherence sqrt(1 - S3^2), S2 = 0
    rng = np.random.default_rng([seed, 1])
    F0, F1 = readout_F
    s3 = 2.0 * (float(np.clip(0.5 + source_bias, 0.02, 0.98)) - 0.5)
    s1 = float(coherence) * math.sqrt(1.0 - s3 * s3)
    out = []
    for s in (s1, 0.0):
        true = (rng.random((n_shots, n_qubits)) > (1.0 + s) / 2.0).astype(int)
        flip = rng.random((n_shots, n_qubits)) < np.where(true == 0, 1.0 - F0, 1.0 - F1)
        out.append(np.where(flip, 1 - true, true))
    return out[0], out[1]
