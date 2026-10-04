# Contributing

Issues and pull requests are welcome.

## Running the tests

```bash
pip install -e ".[dev]"
python -m pytest tests
python tests/run_tests.py    # the same suite with no pytest dependency
```

The NIST known answer tests in `tests/test_nist_kat.py` pin the estimators to NIST's published vectors to
under 1e-9 bits. A change must leave those numbers equal to the reference.

## Conventions

- Every entropy bound is one sided: it is at or below what the source holds. Tests assert
  `attested <= truth`. A change that raises a reported bound needs a proof or a measurement.
- Estimators and tests that implement a standard are ported from its reference implementation and
  validated against published vectors.
- Comments mark headers and functions. The reasons for a design go in the pull request description or the
  documentation.
