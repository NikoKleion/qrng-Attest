# Security

## Scope of reports

Two kinds of report are sent by email, as described below.

- A vulnerability in the code: for example, a file loader that executes content from the file it reads.
- A value that credits more entropy than the source holds: an attested min-entropy, an output length or a
  self-testing bound above the true one, under the assumptions the documentation states. The output of an
  extractor sized by such a value is not uniform, so an error in that direction is treated as a security
  issue.

A value that credits less than the source holds, a crash, or a wrong number in a direction that cannot
lead to output that is not uniform can go in a public issue.

## How to report

Email nikomaxman@gmail.com with the subject line `qrng-attest security`. Include the version
(`pip show qrng-attest`), the input (a captured file, or the Z, X, Y and calibration shots and the readout
matrices), and the command or call that produced the value.

## Supported versions

| Version | Supported |
|---|---|
| 0.1.x | yes |
