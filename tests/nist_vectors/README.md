# NIST test vectors

These files are unmodified copies from NIST's reference implementation of SP 800-90B,
<https://github.com/usnistgov/SP800-90B_EntropyAssessment>.

- `rand1_short.bin`, `rand4_short.bin`, `rand8_short.bin` come from its `bin/` directory.
- `rand1_short.res`, `rand4_short.res`, `rand8_short.res` come from its `cpp/selftest/refdata/` directory and
  hold the reference output that directory publishes for those inputs. Line endings are the only
  difference. `tests/test_nist_kat.py` reads its expected values from them.

NIST-developed software is provided by NIST as a public service and, as a work of the United States
Government, is not subject to copyright protection in the United States (17 U.S.C. 105). It is provided
as is, with no warranty of any kind, express, implied, or statutory, including the implied warranties of
merchantability and fitness for a particular purpose, and NIST assumes no responsibility whatsoever for
its use. These files are not covered by the Apache 2.0 grant that applies to the rest of this
repository. See [NOTICE](../../NOTICE).
