# Tropcc

**Brandon Li**

Finiteness of planar five-body central configurations for generic positive masses in the equal-pair family `(a,a,b,b,c)`, up to similarity.

[Read the paper](paper/tropcc.pdf) · [LaTeX source](paper/tropcc.tex) · [Certificate map](CERTIFICATES.md)

## Verify

Requires Python 3.10 or later. Verification uses exact arithmetic and requires the cddlib command `cddexec_gmp` for the independent fan audit. Singular and gfan are not required. Install cddlib with GMP support before running the complete replay (Homebrew: `brew install cddlib`; Debian/Ubuntu: packages `libcdd-tools python3-z3`). On Linux/aarch64 the requirements file uses the distribution Z3 bindings instead of installing the PyPI source package. The system-site-packages virtual environment below exposes those bindings; the recorded audit used Z3 5.1.0.

```sh
python3 -m venv --system-site-packages .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/fetch_certificates.py
python scripts/verify_all.py
```

The certificate archive is attached to the [v1.0.0 release](https://github.com/brandonlign/tropcc/releases/tag/v1.0.0). The download script checks its SHA-256 digest and every extracted file against `certificate_manifest.json`. An offline copy can be supplied with `--archive /path/to/certificates.tar.xz`. The extracted data occupies approximately 1.3 GB. The full replay includes expensive polynomial elimination.

Successful verification excludes 16,340 of the fan's 16,438 nonzero cones. The remaining 98 cones have 77 distinct rays, all with positive pairing against `(-23,1,1,51,1,1,51,24,13,13)`.

The theorem is relative-generic within the equal-pair mass family. It does not cover every positive mass triple, give a complete explicit exceptional set, or assert finiteness of the entire unaugmented complex distance variety.

## Contents

- `paper/`: manuscript source and PDF.
- `verifier/`: exact certificate, support, symmetry, and coverage checks.
- `data/`: equations, support fan, and separating vectors.
- `companion/`: Laurent proof batches and structural checks.
- `scripts/`: certificate installation and complete replay.

The source and small proof inputs are tracked in Git. Large certificate files are distributed with the release.

## Independent fan audit

The cross-version `fancheck.py` compares two gfan outputs. A separate check uses cddlib with GMP to convert the stored cones to exact half-spaces, then Z3 rational linear arithmetic to ask whether a weight outside their union lies in all 35 tropical hypersurfaces. It checks the exact support and fan symmetries before reducing to six coordinate cases. `UNKNOWN` and timeouts are failures, never completeness certificates.

The complete replay includes this audit. It can also be run separately:

```sh
python verifier/fan_completeness.py --symmetry-break
```

The audit checks an incomplete toy fan as a negative control: removing a ray must produce a witness. Its complete toy fan must return `UNSAT`. Runtime depends on the machine; `--timeout` sets seconds per case.

All six symmetry cases returned `UNSAT` in the checked audit. Input hashes and per-case results are recorded in `data/fan_completeness_audit.json`; `data/fan_completeness.smt2.gz` contains the checked base query. The code rebuilds half-spaces from the original cone rays instead of trusting a cached fan-completeness verdict.
