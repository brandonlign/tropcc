# Tropcc

**Brandon Li**

Finiteness of planar five-body central configurations for generic positive masses in the equal-pair family `(a,a,b,b,c)`, up to similarity.

[Read the paper](paper/tropcc.pdf) · [LaTeX source](paper/tropcc.tex) · [Certificate map](CERTIFICATES.md)

## Verify

Requires Python 3.10 or later. Verification uses exact arithmetic and does not require Singular or gfan.

```sh
python3 -m venv .venv
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
