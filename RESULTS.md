# Results

## Setup

### Synthetic

`synthetic_qrng` generates shots from a device model: a source bias, a readout for each qubit applied as an
asymmetric assignment matrix, optional copying between qubits, and optional drift. Ground truth is known by
construction, so a credited rate can be compared against the entropy the source holds. The
estimators run on the resulting bit stream with no knowledge of the model.

### Hardware

IBM `ibm_marrakesh`, job `dablfc5e36ac739fjmmg`, qubits [0, 1], 2000 shots per circuit, submitted through
qiskit-ibm-runtime with the SamplerV2 primitive. Three circuits: prepare |0> and measure, prepare |1> and
measure, and prepare |+> and measure. The first two calibrate the readout on the same device and in the
same job as the third, which is the source being assessed.

---

## demo

`python -m qrng_attest demo --shots 4000 --qubits 3`, in `results/demo.txt`. A synthetic source with a
known model, showing the certificate format.

```
attested min-entropy = 0.6858 bits/bit   (verdict: PASS)
classical estimate: IID track (source passed the IID checks)
extractable = 9422 uniform bits at epsilon = 2^-64
  0.7852 bits/bit over these 12000 bits taken together, set by the alpha entropy at alpha = 1.229
  long outputs: at least 0.8875 to 0.8893 bits/bit by qubit
source coherence at least 0.970, from X and Y basis shots
```

---

## nist_vectors

`python scripts/nist_vectors.py`, in `results/nist_vectors.txt`. Known answer check against the NIST SP
800-90B reference vectors. The reference column is the output NIST publishes for its `ea_non_iid` on the
same vectors, kept in `tests/nist_vectors/*.res`.

`rand1_short`, binary, per estimator:

| estimator | this package | NIST | difference |
|---|---|---|---|
| most_common_value | 0.9610588257 | 0.9610588257 | 1.1e-16 |
| collision | 0.6914641210 | 0.6914641210 | 4.8e-15 |
| markov | 0.9875961045 | 0.9875961045 | 1.1e-15 |
| compression | 0.6117162048 | 0.6117162048 | 1.1e-16 |
| t_tuple | 0.8676244308 | 0.8676244308 | 0 |
| lrs | 0.9626258038 | 0.9626258038 | 0 |
| multimcw | 0.9526180605 | 0.9526180605 | 0 |
| lag | 0.9433337066 | 0.9433337066 | 0 |
| multimmc | 0.9616166678 | 0.9616166678 | 1.1e-16 |
| lz78y | 0.9614462460 | 0.9614462460 | 0 |

`rand4_short` (alphabet 16) and `rand8_short` (alphabet 256) cover samples wider than one bit, where a source is
assessed as literal symbols and as its bitstring. Every estimator of both tracks, `H_original`,
`H_bitstring` and the assessed value are compared, 20 values per vector:

| vector | H_original | H_bitstring | assessed | worst difference |
|---|---|---|---|---|
| rand4_short | 3.5674726724 | 0.8038720670 | 3.2154882679 | 2.3e-15 |
| rand8_short | 6.6364412871 | 0.7326117181 | 5.8608937445 | 2.3e-11 |

The 2.3e-11 is on the Lag predictor of `rand8_short`, whose estimate is set by the longest run of correct
predictions; the reference evaluates that bound in extended precision. For the four predictors the
reference also prints the number of correct predictions, the longest run and the number of predictions.
These integers are equal on all 20 predictor blocks of the three vectors.

`python scripts/nist_vectors.py <clone> [vector ...]`, in `results/nist_reference_vectors.txt`, runs the
whole assessment on the vectors of a clone of NIST's repository and compares every min-entropy line of
both tracks, `H_original`, `H_bitstring` and the assessed value with `cpp/selftest/refdata/*.res`. All
eleven vectors of the repository at commit 87c104d0, version 1.1.8:

| vector | samples | alphabet | values compared | assessed | NIST | worst difference |
|---|---|---|---|---|---|---|
| biased-random-bits | 1000000 | 2 | 12 | 0.0177665791 | 0.0177665791 | 8.22e-12 |
| biased-random-bytes | 1000000 | 256 | 20 | 0.2577408710 | 0.2577408710 | 1.86e-10 |
| data.pi | 1165666 | 2 | 12 | 0.5695375935 | 0.5695375935 | 4.33e-15 |
| normal | 1000000 | 179 | 20 | 4.1000957831 | 4.1000957831 | 6.39e-14 |
| rand1_short | 10000 | 2 | 12 | 0.6117162048 | 0.6117162048 | 4.77e-15 |
| rand4_short | 10000 | 16 | 20 | 3.2154882679 | 3.2154882679 | 2.33e-15 |
| rand8_short | 10000 | 256 | 20 | 5.8608937445 | 5.8608937445 | 2.27e-11 |
| ringOsc-nist | 1000000 | 2 | 12 | 0.1264457362 | 0.1264457362 | 3.21e-10 |
| truerand_1bit | 1000000 | 2 | 12 | 0.8296770832 | 0.8296770832 | 2.10e-14 |
| truerand_4bit | 1000000 | 16 | 20 | 3.6025063686 | 3.6025063686 | 3.47e-14 |
| truerand_8bit | 1000000 | 256 | 20 | 7.2338614552 | 7.2338614552 | 1.01e-13 |

180 values, worst difference 3.21e-10 bits. On a laptop a vector of a million bits takes 3 to 8 seconds and
one of a million 8 bit samples 68 to 155 seconds.
The reference is NIST's software at one commit; agreement with it is not a NIST validation.

---

## overcredit

`python -m qrng_attest.overcredit`, in `results/overcredit.txt`. A nearly deterministic source, P(0) = 0.980
and true min-entropy 0.0291 bits/bit, read through a characterized symmetric readout of fidelity F, 60000
shots per point.

| F | P(0) observed | passes IID | SP 800-90B estimate | attested | estimate minus 0.0291 |
|---|---|---|---|---|---|
| 0.550 | 0.546 | yes | 0.8580 | 0.0000 | +0.8289 |
| 0.650 | 0.641 | yes | 0.6293 | 0.0170 | +0.6002 |
| 0.750 | 0.740 | yes | 0.4247 | 0.0145 | +0.3956 |
| 0.850 | 0.837 | yes | 0.2498 | 0.0187 | +0.2207 |
| 0.950 | 0.932 | yes | 0.0974 | 0.0247 | +0.0682 |
| 0.999 | 0.979 | yes | 0.0282 | 0.0268 | -0.0010 |

Every observed bit is an independent draw at bias `p_obs = F p + (1 - F)(1 - p)`, so the stream is
IID and passes the IID screen at every point. The SP 800-90B estimate is 0.858 bits/bit at
F = 0.55. The attested value is at or below 0.0291 in every row.

The file also carries the crosstalk and miscalibration experiments, and the coherence experiment: a source
whose Z outcomes are uniform at every coherence c, 60000 shots in each of the Z, X and Y bases, readout
fidelity 0.99 calibrated on as many shots. The true min-entropy is eq. 7 of Fiorentino et al.

| coherence | true min-entropy | SP 800-90B estimate | population value | coherence bound |
|---|---|---|---|---|
| 1.000 | 1.0000 | 0.9834 | 0.9801 | 0.7277 |
| 0.999 | 0.9369 | 0.9834 | 0.9801 | 0.7220 |
| 0.990 | 0.8096 | 0.9834 | 0.9801 | 0.6740 |
| 0.950 | 0.6080 | 0.9834 | 0.9801 | 0.5410 |
| 0.900 | 0.4781 | 0.9834 | 0.9801 | 0.4372 |
| 0.500 | 0.1000 | 0.9834 | 0.9801 | 0.0887 |
| 0.000 | 0.0000 | 0.9834 | 0.9801 | 0.0000 |

---

## hardware_ibm_marrakesh

Measured, in `results/hardware_ibm_marrakesh.txt` with raw values in the accompanying `.json`.

Readout calibration, from the |0> and |1> preparations in the same job:

| qubit | P(read 1 \| prep 0) | P(read 1 \| prep 1) | F0 | F1 |
|---|---|---|---|---|
| 0 | 0.0050 | 0.9945 | 0.9950 | 0.9945 |
| 1 | 0.0045 | 0.9900 | 0.9955 | 0.9900 |

The |+> circuit gives P(1) = 0.5115 and 0.5135 on the two qubits. The certificate:

```
attested min-entropy = 0.8607 bits/bit   (verdict: FLAGGED)
classical estimate: IID track (source passed the IID checks)
extractable = 3316 uniform bits at epsilon = 2^-64
action: q1: readout error sets the attested rate (0.86 bits/bit); recalibrate the discriminator or exclude this qubit
```

The verdict is FLAGGED because the readout error of qubit 1 sets the attested rate. The job holds Z basis
shots only, so 0.8607 is the population value and holds for a pure source state.

Coverage: one backend, two qubits, 2000 shots, one calibration cycle.

Shots taken together, a known coherence and no bias, bits per shot certified at a smoothing of 2^-65:

| coherence | one shot | n = 20000 | n = 1000000 | n = 100000000 | long outputs |
|---|---|---|---|---|---|
| 1.000 | 1.0000 | 0.9934 | 0.9999 | 1.0000 | 1.0000 |
| 0.998 | 0.9116 | 0.9523 | 0.9842 | 0.9882 | 0.9886 |
| 0.995 | 0.8627 | 0.9287 | 0.9688 | 0.9742 | 0.9748 |
| 0.990 | 0.8096 | 0.8994 | 0.9472 | 0.9539 | 0.9546 |
| 0.950 | 0.6080 | 0.7520 | 0.8202 | 0.8302 | 0.8313 |
| 0.900 | 0.4781 | 0.6264 | 0.7011 | 0.7124 | 0.7136 |
| 0.500 | 0.1000 | 0.1287 | 0.1796 | 0.1878 | 0.1887 |

One shot is eq. 7. The columns for n shots are Theorem 7 of Tomamichel, Colbeck and Renner with the alpha
entropy of one shot at the best alpha in (1, 2]. Long outputs is the conditional von Neumann entropy. For a
pure state the output is n coin flips, whose smooth min-entropy has an exact value; the bound is below it
by 0.002 to 0.02 bits per bit at n of 1000 and 5000.

---

## hardware_coherence_ibm_fez

Measured, by `scripts/coherence_hardware.py`, in `results/hardware_coherence_ibm_fez.txt`, with the
outcomes of every shot in `_shots.npz`, the job record in `.json` and the executed circuits in `_circuits.qpy`. IBM
`ibm_fez`, job `db0rsb3id5ic73eq8em0`, qubits 6, 124, 144 and 50, 20000 shots per circuit, 28 quantum seconds
on the open plan. Five circuits in one job: prepare |0> and measure, prepare |1> and measure, and |+> measured
in Z, in X and in Y. Dynamical decoupling and twirling off. The qubits are the four of least readout plus sx
error with no two coupled.

| qubit | F0 | F1 | P(0) in Z | coherence | at least | population value | coherence bound |
|---|---|---|---|---|---|---|---|
| 6 | 0.9986 | 0.9919 | 0.5076 | 0.9982 | 0.9930 | 0.9582 | 0.8388 |
| 124 | 0.9997 | 0.9929 | 0.4995 | 0.9983 | 0.9956 | 0.9605 | 0.8702 |
| 144 | 0.9996 | 0.9958 | 0.5004 | 0.9982 | 0.9957 | 0.9681 | 0.8721 |
| 50 | 0.9956 | 0.9935 | 0.4947 | 0.9965 | 0.9843 | 0.9518 | 0.7653 |

The certificate:

```
attested min-entropy = 0.7653 bits/bit   (verdict: PASS)
classical estimate: IID track (source passed the IID checks)
extractable = 75020 uniform bits at epsilon = 2^-64
  0.9377 bits/bit over these 80000 bits taken together, set by the alpha entropy at alpha = 1.129
  long outputs: at least 0.9337 to 0.9778 bits/bit by qubit
source coherence at least 0.984, from X and Y basis shots
```

The measured coherence is 0.9965 to 0.9983, and the attested value is the coherence bound of qubit 50 at its
lower limit, 0.9843. That limit takes the error seen in the calibration of |0>, 0.0044 on this
qubit, as a possible fault of the prepared state. Eq. 7 at the measured coherences themselves gives 0.8845
to 0.9175. The population value, 0.95 to 0.97, is above the bound on every qubit.

The attested value is for one shot and gives 61096 bits by the leftover hash lemma. The 80000 bits taken
together have a smooth min-entropy of at least 75148 bits at smoothing 2^-65 (alpha = 1.129), which
gives the 75020 of the certificate. Long outputs at these lower limits on the coherence approach
0.9337 to 0.9778 bits per bit, by qubit.

Coverage: one backend, four qubits, one job, one calibration cycle.

---

## hardware_coherence_rounds_ibm_fez

Measured, by `scripts/coherence_rounds.py`, in `results/hardware_coherence_rounds_ibm_fez.txt`, with the
outcomes of every shot of every circuit in `_shots.npz`, the job record in `.json` and the executed circuits in
`_circuits.qpy`. IBM `ibm_fez`, job `db1693rid5ic73eqnak0`, qubits 2, 130, 124, 86, 144, 14, 64, 23, 110, 90, 75, 4, 111
quantum seconds on the open plan against an estimate of 106 and a cap of 169. One job of five rounds;
a round is prepare |0>, prepare |1>, and |+> measured in Z, in X and in Y, at 20000 shots a circuit and 4000
for Z. Dynamical decoupling and twirling off. The qubits are the twelve of least readout plus sx error with
no two coupled.

All rounds pooled:

| qubit | F0 | F1 | P(0) in Z | coherence | at least | population value | coherence bound | eq. 7 at the measured coherence |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.9987 | 0.9944 | 0.5062 | 0.9996 | 0.9964 | 0.9611 | 0.8821 | 0.9609 |
| 130 | 0.9991 | 0.9897 | 0.4973 | 0.9995 | 0.9969 | 0.9511 | 0.8899 | 0.9539 |
| 124 | 0.9991 | 0.9861 | 0.4981 | 0.9990 | 0.9964 | 0.9480 | 0.8826 | 0.9383 |
| 86 | 0.9987 | 0.9926 | 0.5012 | 0.9988 | 0.9953 | 0.9672 | 0.8666 | 0.9310 |
| 144 | 0.9998 | 0.9918 | 0.4928 | 0.9991 | 0.9981 | 0.9413 | 0.9133 | 0.9410 |
| 14 | 0.9980 | 0.9935 | 0.5001 | 0.9977 | 0.9927 | 0.9661 | 0.8360 | 0.9047 |
| 64 | 0.9991 | 0.9934 | 0.4935 | 0.9987 | 0.9962 | 0.9459 | 0.8793 | 0.9276 |
| 23 | 0.9990 | 0.9935 | 0.4982 | 0.9972 | 0.9944 | 0.9595 | 0.8556 | 0.8963 |
| 110 | 0.9976 | 0.9929 | 0.4979 | 0.9987 | 0.9928 | 0.9596 | 0.8369 | 0.9285 |
| 90 | 0.9994 | 0.9921 | 0.4758 | 0.9994 | 0.9975 | 0.8952 | 0.9016 | 0.9504 |
| 75 | 0.9987 | 0.9879 | 0.4918 | 0.9987 | 0.9952 | 0.9333 | 0.8648 | 0.9284 |
| 4 | 0.9969 | 0.9538 | 0.5071 | 0.9982 | 0.9906 | 0.9269 | 0.8149 | 0.9149 |

The coherence measured in each round, from that round's own calibration:

| qubit | 1 | 2 | 3 | 4 | 5 | spread |
|---|---|---|---|---|---|---|
| 2 | 1.0000 | 0.9999 | 0.9986 | 0.9998 | 0.9994 | 0.0014 |
| 130 | 0.9985 | 0.9996 | 0.9984 | 1.0000 | 1.0000 | 0.0016 |
| 124 | 0.9987 | 0.9989 | 0.9993 | 1.0000 | 0.9984 | 0.0016 |
| 86 | 0.9982 | 0.9986 | 0.9998 | 0.9982 | 0.9993 | 0.0016 |
| 144 | 0.9993 | 0.9997 | 0.9989 | 0.9991 | 0.9987 | 0.0010 |
| 14 | 0.9976 | 0.9977 | 0.9987 | 0.9968 | 0.9977 | 0.0020 |
| 64 | 0.9979 | 0.9986 | 0.9984 | 0.9997 | 0.9990 | 0.0018 |
| 23 | 0.9980 | 0.9976 | 0.9973 | 0.9970 | 0.9964 | 0.0017 |
| 110 | 0.9980 | 0.9993 | 0.9992 | 0.9994 | 0.9978 | 0.0016 |
| 90 | 0.9986 | 0.9998 | 0.9998 | 0.9998 | 0.9990 | 0.0012 |
| 75 | 0.9979 | 0.9996 | 0.9994 | 0.9985 | 0.9982 | 0.0017 |
| 4 | 0.9980 | 0.9979 | 0.9982 | 0.9997 | 0.9971 | 0.0026 |

The coherence bound after the first rounds, pooled:

| qubit | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| 2 | 0.8794 | 0.8833 | 0.8800 | 0.8814 | 0.8821 |
| 130 | 0.8699 | 0.8849 | 0.8827 | 0.8879 | 0.8899 |
| 124 | 0.8566 | 0.8722 | 0.8768 | 0.8846 | 0.8826 |
| 86 | 0.8442 | 0.8572 | 0.8648 | 0.8638 | 0.8666 |
| 144 | 0.8940 | 0.9121 | 0.9104 | 0.9131 | 0.9133 |
| 14 | 0.8232 | 0.8304 | 0.8344 | 0.8358 | 0.8360 |
| 64 | 0.8527 | 0.8654 | 0.8711 | 0.8754 | 0.8793 |
| 23 | 0.8456 | 0.8544 | 0.8558 | 0.8554 | 0.8556 |
| 110 | 0.8071 | 0.8270 | 0.8320 | 0.8354 | 0.8369 |
| 90 | 0.8755 | 0.8921 | 0.8994 | 0.9034 | 0.9016 |
| 75 | 0.8411 | 0.8590 | 0.8620 | 0.8638 | 0.8648 |
| 4 | 0.7939 | 0.8009 | 0.8083 | 0.8128 | 0.8149 |

The certificate:

```
attested min-entropy = 0.8149 bits/bit   (verdict: PASS)
classical estimate: IID track (source passed the IID checks)
extractable = 229744 uniform bits at epsilon = 2^-64
  0.9573 bits/bit over these 240000 bits taken together, set by the SP 800-90B values
  long outputs: at least 0.9568 to 0.9890 bits/bit by qubit
source coherence at least 0.991, from X and Y basis shots
```

The bound of a qubit is 0.8149 to 0.9133 at 100000 shots for the calibration and for X and Y. The measured
coherences are 0.9972 to 0.9996, and eq. 7 at those values is 0.8963 to 0.9609. The spread of a qubit's
coherence over the five rounds is at most 0.0026, of the size of the sampling error of one round, so the
job shows no drift between rounds.

The output length is set by the SP 800-90B values of the qubits, 0.9070 to 0.9737 at 20000 Z shots each, which
sum to 0.9578 per bit; the alpha entropy of the shots taken together allows more. Qubit 90 reads 0 in
47.58 percent of its Z shots and has the least of those values. Long outputs at the lower limits on the
coherence approach 0.9568 to 0.9890 bits per bit, by qubit.

Coverage: one backend, twelve qubits, one job, one calibration cycle.
