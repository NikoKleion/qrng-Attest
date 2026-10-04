# Using qrng-attest

## Overview

qrng-attest estimates the min-entropy of a quantum random source and extracts uniform bits at that rate.
The observed output includes measurement noise, correlation between qubits and device drift, which the
attestation checks separately.

## Checks

The entropy math is checked against NIST's published test vectors and a test suite over synthetic data
and captured files. Reading a bit stream runs the same path whether the bits come from a file or a device.

The backend path runs on IBM hardware through qiskit-ibm-runtime: job submission, the bit order of each qubit,
readout calibration, result parsing, and the certificate. It has run on ibm_fez with four and with
twelve qubits in the Z, X and Y bases and on ibm_marrakesh with two qubits in the Z basis, so on another
backend check the bit order of each qubit and the calibration fields that backend exposes.

## Terms

- **Min-entropy.** H_min = -log2 max_x P(x), in bits per sample (NIST SP 800-90B, section 3.1). A fair coin
  gives 1 bit per flip. This is the quantity the package reports.
- **SP 800-90B assessment.** The non-IID estimators of section 6.3 and the IID track of section 5 of NIST
  SP 800-90B. IID means independent and identically distributed samples. The assessed value is the minimum
  over the estimators; a source that passes the IID track is assessed by the most common value estimate
  (section 6.1).
- **Readout error, assignment matrix.** M[a, b] = P(read a | prepared b) for one qubit, a 2x2 matrix
  measured from |0> and |1> preparations. Readout error adds randomness that the quantum state did not
  produce. The same matrix is the calibration matrix of measurement error mitigation (Bravyi et al., Phys.
  Rev. A 103, 042605 (2021)).
- **Attested min-entropy.** The smallest of three values, each at a 99 percent one sided confidence
  margin: the SP 800-90B assessment of the bits; the min-entropy of the Z populations after the assignment
  matrix is inverted out of the observed frequencies, which is the min-entropy of a pure source state
  (Fiorentino et al., Phys. Rev. A 75, 032334 (2007), eq. 4); and the coherence bound. The checks of
  correlation between qubits and of drift follow.
- **Coherence bound.** For a source state with coherence c = sqrt(S1^2 + S2^2), where S1 and S2 are its X and
  Y expectation values, the min-entropy of the Z outcome against an adversary who holds a purification is
  -log2((1 + sqrt(1 - c^2)) / 2) (Fiorentino et al., eq. 7). The package measures S1 from shots taken after a
  Hadamard and S2 from shots taken after an S-dagger and a Hadamard, corrects both for readout with each
  fidelity taken between its measured value less a margin and 1 (an error seen in the calibration may be a
  fault of the prepared state), and uses a lower bound on c. A state with no coherence, such as a classical
  mixture of |0> and |1>, gets 0 whatever its Z statistics. Without X and Y shots the certificate carries a
  note that its value assumes a pure source state.
- **Extraction.** A randomness extractor maps n input bits with min-entropy k to ell output bits within
  epsilon of uniform. The package uses Toeplitz hashing (Krawczyk, CRYPTO 1994) sized by the leftover hash
  lemma, ell = k - 2 log2(1 / epsilon) + 2 (Impagliazzo, Levin and Luby, STOC 1989).
- **Shots taken together.** With X and Y shots, the output length is the larger of that rule at the attested
  one shot value and ell = K - 2 log2(1 / epsilon), where K is a lower bound on the smooth min-entropy of the
  whole output at smoothing epsilon / 2: the largest over alpha in (1, 2] of n times the summed alpha
  entropies of the qubits less log2(8 / epsilon^2) / (alpha - 1) (Tomamichel, Colbeck and Renner,
  arXiv:0811.1221, Theorem 7; Tomamichel et al., arXiv:1002.2436, Theorem 6). K is capped at n times the sum
  of the qubits' SP 800-90B values, and at the bound for correlated qubits. `certificate.plan.order` is the
  alpha used, or None where the one shot rule gives more; `attestation.per_qubit[i]["h_long_run"]` is the rate
  long outputs approach. The shots are taken as identical and independent, as for the one shot value.
- **Self-testing bound.** `self_testing_entropy(ones, shots)` takes the counts of outcome 1 for four
  preparations x and two measurements y, and returns the dimension witness W, the absolute value of the
  determinant with rows (p(1|0,y) - p(1|1,y), p(1|2,y) - p(1|3,y)) for y = 0, 1, a lower limit on W over the
  Wilson intervals of the eight proportions, and -log2 of the largest guessing probability at that limit, (1 +
  sqrt((1 + sqrt(1 - W^2)) / 2)) / 2 (Lunghi et al., arXiv:1410.2790, eqs. 4 and 5). The entropy is per round,
  averaged over the eight settings taken equally often. It assumes a qubit, preparation and measurement
  devices that are independent, identical rounds, and devices that do not use the choice of setting; it does
  not use the readout calibration, or the gates that change the basis being as specified.
  `qiskit_bridge.selftest_circuits(nq)` gives the eight circuits: |0>, |1>, |+> and |-> measured in Z and in
  X, for which W = 1 and the bound is 0.2284.
- **Certificate, PASS, FLAGGED.** The printed report. PASS needs X and Y shots, an attested value above
  0.01, and every check passing. FLAGGED means one of these failed; an `action` line names it. `advice`
  lines say what sets the rate and are not failures.
- **Purification, pure source state.** A pure state is one that is not mixed with anything; a mixed state can
  be written as part of a larger pure state, and the rest of that larger state is the purification. The
  coherence bound assumes an adversary holds it, which is the most an adversary could know about the qubit.
- **Alpha entropy, smooth min-entropy, conditional von Neumann entropy.** The alpha entropy is the Renyi
  entropy of order alpha of the outcome given the adversary's system; at alpha near 1 it is the conditional
  von Neumann entropy, the rate long outputs approach. The smooth min-entropy at smoothing epsilon is the
  largest min-entropy of any state within epsilon of the real one; it is what an extractor's output length
  depends on.
- **Health tests, drift test, readout check.** The health tests are the repetition count and adaptive
  proportion tests of SP 800-90B section 4.4 on each qubit's bits. The drift test compares P(0) across
  twenty windows of the first qubit's stream. The readout check fails a qubit whose value corrected for readout
  is more than 0.1 bits under its observed one.
- **Literal and bitstring, `H_original` and `H_bitstring`.** Samples wider than one bit are assessed twice, as
  symbols (literal) and as the string of their bits; the result is the least of the symbol width,
  `H_original` and the width times `H_bitstring` (SP 800-90B section 3.1.3). Predictors are four of the
  estimators of section 6.3: MultiMCW, Lag, MultiMMC and LZ78Y.
- **Conditioning, restart tests, CAVP and CMVP.** Conditioning is any processing after the raw samples; the
  package assesses the samples before it. The restart tests of SP 800-90B section 3.1.4 are not included.
  CAVP and CMVP are NIST's validation programs; the package is not part of one.
- **Recorded device model, open plan, quantum seconds.** A fake backend of qiskit-ibm-runtime is a snapshot of
  a device's calibration that runs locally. The open plan is IBM Quantum's free plan; quantum seconds are the
  device time a job is charged, which the plan limits per period.

## Install

```bash
pip install git+https://github.com/NikoKleion/qrng-attest
```

Not on PyPI. Install from the repository, or from a clone with `pip install -e .`.

Requires Python 3.11 or newer and numpy. Qiskit is optional and only needed to certify a backend. Extras
install it:

```bash
pip install -e ".[qiskit]"    # backend bridge
pip install -e ".[runtime]"   # run on IBM hardware through qiskit-ibm-runtime
pip install -e ".[aer]"       # run under a local noise model
pip install -e ".[figures]"   # redraw the figures in README.md
pip install -e ".[sp80022]"   # NIST SP 800-22 battery through nistrng
```

Assessing a captured file or an array needs only numpy.

## Demo

This runs a simulated source and prints a certificate. It needs no hardware.

```bash
python -m qrng_attest demo --qubits 4
```

The synthetic source is a pure |+> state read through a readout of fidelity 0.99 and 0.98, so its true
min-entropy is 1 bit per bit. The certificate's value is what the shots and that readout let the package
certify.

Flags add faults to the synthetic source:

```bash
python -m qrng_attest demo --qubits 4 --crosstalk 0.5 --drift 0.2
```

`--coherence 0.9` sets the coherence of the synthetic source state, and `--no-tomography` leaves out its X
and Y basis shots.

## Assess a file of captured bits

A file of raw output is assessed directly.

```bash
python -m qrng_attest capture mydata.bin
```

The file is read as a bit stream by default. If each byte is a symbol, add `--bytes` (`mode="bytes"` in
Python). The report gives the min-entropy, which estimator set the limit, and how many uniform bits it
gives. The number to quote is the `min-entropy =` line, `report.min_entropy`.

SP 800-90B asks for 1 000 000 samples. On a shorter stream the non-IID estimators report less than the
source holds: 20000 fair bits assess at about 0.64 to 0.85 bits per bit without the IID track. The
chi-square independence test of the IID track holds its one in a thousand rejection rate on 8 bit samples
at a million samples; on uniform bytes it rejects 65 percent of streams of 20000 samples, 8 percent at
50000 and 3 percent at 100000, and the report then gives the non-IID value. On bits it holds its rate at
every size tried. The IID track is not compared with NIST's `ea_iid`.

From Python:

```python
from qrng_attest import assess

report = assess("mydata.bin", mode="bits", run_iid=True)
print(report)
print(report.min_entropy)
```

## From arrays of shots to uniform bits

Every array of shots has shape (shots, qubits) and holds 0 and 1. A source takes the Z shots, which become
the output; X and Y shots of the same preparation, which may be more or fewer than the Z shots; and the
readout, as one 2x2 assignment matrix for every qubit or a list of one per qubit, with the number of shots
it was measured from.

```python
from qrng_attest import (QuantumArraySource, assignment_matrix, calibrate_readout, certify_source,
                         extract_from_source, synthetic_qrng, synthetic_tomography)

z, matrices = synthetic_qrng(8000, n_qubits=3, readout_F=(0.99, 0.98), seed=1)         # Z shots, one matrix a qubit
x, y = synthetic_tomography(32000, n_qubits=3, coherence=0.995, readout_F=(0.99, 0.98), seed=1)

source = QuantumArraySource(z, readout=matrices, calibration_shots=32000, x_shots=x, y_shots=y)
cert = certify_source(source)
print(cert)
bits, cert, seed = extract_from_source(source)       # len(bits) == cert.plan.output_bits
```

- `synthetic_qrng(n_shots, n_qubits=4, source_bias=0.0, readout_F=(0.99, 0.98), crosstalk=0.0, drift=0.0,
  seed=0)` returns the Z shots of a |+> source and the list of assignment matrices. `source_bias` moves P(0)
  off one half, `crosstalk` is the probability that a qubit copies its neighbour, `drift` is the change of
  the bias over the run.
- `synthetic_tomography(n_shots, n_qubits=4, source_bias=0.0, coherence=1.0, readout_F=(0.99, 0.98), seed=0)`
  returns the X shots and the Y shots of the same state at the given coherence.
- With measured data, `calibrate_readout(prep0_outcomes, prep1_outcomes)` returns the assignment matrix of one
  qubit and the number of calibration shots, and `assignment_matrix(F0, F1)` builds one from two fidelities.
- A source that fails: `crosstalk=0.5` in `synthetic_qrng` is FLAGGED for correlated qubits, and
  `coherence=0.0` in `synthetic_tomography` is attested at zero. `extract_from_source` raises a RuntimeError
  for a FLAGGED source; `require_pass=False` extracts anyway.

## Certify a Qiskit backend

This prepares |+> on each qubit of the backend and measures it in the Z, X and Y bases, measures the
readout error on the same device, attests the result, and reports how many uniform bits it supports.
`tomography=False` leaves out the X and Y circuits, and the certificate then assumes a pure source state.
It runs the IID track by default; on a source close to uniform that passes it, the coherence bound sets the
attested rate. Pass `run_iid=False` to use the non-IID estimate.

```python
from qiskit_ibm_runtime.fake_provider import FakeManilaV2
from qrng_attest import certify_backend

cert = certify_backend(FakeManilaV2(), shots=8192, simulate=True)
print(cert)
print(cert.rate)            # attested bits per bit
print(cert.passed)          # True if the source cleared every check
print(cert.key_cost(256))   # raw bits needed for a 256 bit key
```

`qubits` is a list of physical qubits; left out, the first five of the device are read. With
`simulate=True` the circuits run under a local noise model built from the backend, and the readout matrix
is taken from that model, so no calibration circuits run and the certificate notes that the matrix is
taken as exact. With `simulate=False` the circuits are run by the backend itself, with the calibration
circuits: a fake backend runs them locally and is labelled synthetic; a backend object of
`QiskitRuntimeService` is a live device, which needs the runtime extra and waits in the queue.

On a live device `certify_backend` submits at once: it does not look at the account's plan or allowance and
sets no time cap. Its five circuits of `shots` shots and two calibration circuits take about a quarter of
a millisecond a shot. The scripts of the next section are the path with those checks.

From the command line, `qrng-attest backend <name>` takes the name of a fake backend, for example
`manila`, `fez` or `brisbane`, and `--qubits 0,1,2` chooses the qubits. It does not reach a live device.

## Run on a live IBM device

`certify_backend(backend, shots=...)` with a backend object from `QiskitRuntimeService` submits the circuits
and waits for them. The two scripts below take the same circuits in three separate steps, so that the plan
and its cost are seen before anything is submitted. Both need the `runtime` extra and a saved account;
`python scripts/setup_account.py ibm_creds.txt` saves one from a local file that holds a line `token=...`,
so the token is not typed on a command line.

```bash
python scripts/coherence_hardware.py dry ibm_fez     # 4 qubits, 5 circuits of 20000 shots; nothing is submitted
python scripts/coherence_hardware.py submit
python scripts/coherence_hardware.py collect

python scripts/coherence_rounds.py dry ibm_fez       # 12 qubits, 5 rounds: 25 circuits, 420000 shots
python scripts/coherence_rounds.py submit            # prints the plan and any reason to stop
python scripts/coherence_rounds.py submit go         # submits, then watches the job to its end
python scripts/coherence_rounds.py collect
```

- `dry` reads the account's plan and allowance, picks the qubits of least readout plus sx error with no two
  coupled, transpiles the circuits, checks that each is read on those qubits and that the basis changes are
  still in them, estimates the execution time, rehearses the circuits under the backend's noise model in
  Aer, and writes the plan under `build/`. It submits nothing.
- Device time is about the number of shots times the circuit duration plus the repetition delay, 0.25 ms a
  shot on ibm_fez, whatever the number of qubits. The two jobs in [RESULTS.md](RESULTS.md) used 28 and 111
  seconds.
- `submit` stops unless the account has exactly one instance, on the open plan with free pricing, with at
  least twice the estimate left of its allowance, and it sets the job's `max_execution_time`. In
  `coherence_rounds.py` it also stops for a plan older than an hour, a plan that is not the script's own
  job, an estimate above 117 s, or an allowance that would fall under 120 s if the job ran to its cap; it
  submits only on the word `go`, and then cancels the job if it runs 20 s past its cap.
- `collect` saves the outcome of every shot and computes the certificate from the saved file, so a result
  is recomputed without the account.

A plan goes stale: the qubits of least error change with the device's calibration. Run `dry` again before
`submit` when time has passed.

## Self-testing bound

```python
from qiskit_ibm_runtime.fake_provider import FakeManilaV2
from qrng_attest import self_testing_entropy
from qrng_attest.qiskit_bridge import selftest_circuits, selftest_on_backend

circuits = selftest_circuits(2)       # a dict keyed (x, y): preparation x in 0..3 (|0>, |1>, |+>, |->), measurement y in 0..1 (Z, X)
results = selftest_on_backend(FakeManilaV2(), shots=4000, qubits=[0, 1])   # the eight circuits under the device model
print(results[0]["witness"], results[0]["witness_low"], results[0]["h_min"])

N = 20000                             # or from counts: ones[x][y] outcomes 1 in shots[x][y] rounds, both 4 by 2
ones = [[0, N // 2], [N, N // 2], [N // 2, 0], [N // 2, N]]
print(self_testing_entropy(ones, [[N, N]] * 4)["h_min"])
```

`self_testing_entropy` returns a dict: `witness`, the dimension witness W of the observed proportions;
`witness_low`, its lower limit; `h_min`, the min-entropy per round at that limit; `p`, the proportions. A round
is one preparation and one measurement, and the entropy is an average over the eight settings taken equally
often.

## Readout calibration

The package measures the readout error by preparing known 0 and 1 states and counting how often they read
back wrong. `BackendV2.target` does not carry asymmetric readout fidelities, and the per direction numbers
that `backend.properties()` does carry are from the vendor's own calibration cycle, which can be older
than the run. Measuring in the same job gives the matrix the provenance of the source circuit. The
calibration shot count sets the uncertainty of the measured matrix, and the reported entropy is the worst
case over that uncertainty. More calibration shots give a smaller uncertainty.

The calibration cannot tell a readout error from a fault of the prepared state, such as a qubit that was
not reset. For the coherence bound each fidelity is therefore taken anywhere from its measured value less
its margin up to 1, and the least resulting coherence is used.

`certify_backend` runs this step. On its own:

```python
from qrng_attest.qiskit_bridge import self_calibrate_backend

matrices, shots = self_calibrate_backend(backend, shots=4000)
```

## Qubit selection

One bad qubit lowers the whole device's rate. The package reports which subset of qubits gives the most
total bits per shot.

```python
from qrng_attest import qubit_subset

result = qubit_subset(source)
print(result["best_qubits"])          # the subset to keep
print(result["dropped"])              # the qubits to leave out
print(result["total_bits_per_shot"])  # throughput with that subset
```

## Extraction

Once a source is certified, extract uniform bits from it.

```python
from qrng_attest import extract_from_source

bits, cert, seed = extract_from_source(source)
```

`bits` is a numpy array of 0s and 1s close to uniform. `seed` is the public seed the extractor used; it is
not secret and can be reused. To size the run before extracting, use the planner:

```python
from qrng_attest import extraction_plan

plan = extraction_plan(n_input=1000000, h_per_symbol=0.68)
print(plan)
print(plan.bits_for(256))   # input samples needed for a 256 bit key
```

## SP 800-22 battery

With the `sp80022` extra installed, the fifteen tests of NIST SP 800-22 rev 1a run on a 0/1 array through
the nistrng package.

```python
from qrng_attest import sp800_22

rows = sp800_22.run_battery(bits)
print(sp800_22.report(rows))
```

A test the stream is too short for is reported as not eligible. `names` runs a subset, for example
`run_battery(bits, names=["monobit", "runs"])`.

## Certificate output

```
QRNG certificate: ibm_fez QRNG  [measured]
  20000 shots x 12 qubits
  attested min-entropy = 0.8149 bits/bit   (verdict: PASS)
  classical estimate: IID track (source passed the IID checks)
  extractable = 229744 uniform bits at epsilon = 2^-64
    0.9573 bits/bit over these 240000 bits taken together, set by the SP 800-90B values
    long outputs: at least 0.9568 to 0.9890 bits/bit by qubit
  a 256-bit key costs 469 raw bits
  source coherence at least 0.991, from X and Y basis shots
```

- The first line names the source and, in brackets, where its shots come from: `[measured]` on a device,
  `[synthetic]` from a model or a generator, `[external]` from a file. The second gives the Z shots per
  qubit and the number of qubits. The raw bits are their product.
- `attested min-entropy` is a lower bound on the min-entropy of one raw bit, the least over the qubits. Its
  verdict is PASS when the health tests, the independence of the qubits, the drift test and the readout
  check all pass. FLAGGED means one failed; an `action` line follows for each, naming the qubit and a remedy:
  recalibrate, drop a correlated qubit, shorten the time between calibrations, or take more X, Y and
  calibration shots.
- `classical estimate` says which SP 800-90B track gave each qubit's classical value: the IID track for a
  qubit that passes the IID checks, the non-IID estimators otherwise.
- `extractable` is the number of uniform bits the raw bits of this run give at the stated epsilon, and the
  length `extract_from_source` returns. The line under it gives that length per raw bit and the term that
  sets it: the alpha entropy of the shots taken together with its alpha, the SP 800-90B values, or the bound
  for correlated qubits. `long outputs` is the rate a longer run of each qubit approaches. These two lines
  are present when X and Y shots are.
- `a 256-bit key costs` is the number of raw bits that give 256 uniform bits at the attested one shot value.
- `source coherence at least` is the least lower limit on the coherence over the qubits. Without X and Y
  shots the line is absent, the verdict is FLAGGED, and the value is an upper bound unless the source
  state is pure.
- `best qubit subset` appears when leaving out some qubits gives more bits per shot, with the qubits to
  drop and the bits per shot before and after.
- `advice` lines say what sets the attested value and how to raise it. `note` lines state what the run
  takes as given: a synthetic source, a readout matrix taken as exact, missing X and Y shots.

## Security parameter epsilon

Extraction takes a security parameter epsilon, the largest allowed distance between the output and a
uniform stream. The default is 2 to the power of minus 64. A smaller epsilon gives output closer to
uniform and fewer bits. Pass it to the extract and plan functions.

```python
plan = extraction_plan(n_input=1000000, h_per_symbol=0.68, epsilon=2.0 ** -32)
```

## IID permutation budget

The `assess` IID check has two modes. The default is a reduced shuffle budget, which differs from the full
one only by denying the IID credit more often. For the reference budget of 10000 shuffles, pass `--iid-full` on the command line
or `iid_full=True` in code.

The certificate path (`certify_source`, `certify_backend`) runs the full 10000-shuffle budget for its
IID check of each qubit. The reduced budget uses a fixed pass threshold calibrated for the full budget, so at
a smaller budget it fails clean data and denies the credit. The permutation test stops early on a clean
stream, so the full budget is cheap on good data. A stream that fails the chi-square or LRS check, which run first,
skips the permutation test, so a bad source stays fast.

## Limitations

It assumes a trusted and characterized measurement, the gates that change the basis included, and the
same source state in every shot. Fiorentino et al. state the same two limits for their method: the state
must not change between the tomography and the bits, and the security covers individual attacks. The
output length for shots taken together assumes the shots are identical and independent, with an
adversary whose system is a product over the shots. The package reimplements the NIST SP 800-90B tests
and equals NIST's reference output on its vectors. It is not device-independent and is not a validation
under CAVP or CMVP.

Each bound is a 99 percent one sided bound on its own quantity; the margins are not combined into one
confidence level, and they are not part of epsilon. The self-testing bound assumes a qubit, preparation
and measurement devices that are independent, identical rounds, and devices that do not use the choice
of setting; it has run under backend noise models only. The three bases have run on one backend.
