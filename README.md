# Tropcc

**Brandon Li**

Finiteness of planar five-body central configurations for generic positive masses in the equal-pair family `(a,a,b,b,c)`, up to similarity.

[Read the paper](paper/tropcc.pdf) · [LaTeX source](paper/tropcc.tex) · [Certificate map](CERTIFICATES.md)

## Verify

Requires Python 3.10 or later, and `yices-smt2` from Yices 2.7 for the second solver audit (Homebrew: `brew install yices2`; other platforms: the official Yices distribution). Verification uses exact arithmetic and requires the cddlib command `cddexec_gmp` for the independent fan audit. Singular and gfan are not required. Install cddlib with GMP support before running the complete replay (Homebrew: `brew install cddlib`; Debian/Ubuntu: packages `libcdd-tools python3-z3`). On Linux/aarch64 the requirements file uses the distribution Z3 bindings instead of installing the PyPI source package. The system-site-packages virtual environment below exposes those bindings; the recorded audit used Z3 5.1.0.

```sh
python3 -m venv --system-site-packages .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/fetch_certificates.py
python scripts/verify_all.py
```

The certificate archive is attached to the [v1.2.3 release](https://github.com/brandonlign/tropcc/releases/tag/v1.2.3). The download script checks its SHA-256 digest and every extracted file against `certificate_manifest.json`. An offline copy can be supplied with `--archive /path/to/certificates.tar.xz`. The extracted data occupies approximately 1.3 GB. The full replay includes expensive polynomial elimination.

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

## Fan structure and mass specialization

`python verifier/fan_structure.py` independently reconstructs all 16,438 cones as Newton normal regions, enumerates their facets in exact rational arithmetic, and checks closure under faces, primitive extreme rays, zero lineality, and invariance of the full fan and all supports. It checks 44,998 facet incidences without gfan, cddlib, or SMT. This verifies the face relations used by the modular-star argument.

`data/cone_accounting.json` maps every covered cell to its proof class and source representative, with all symmetry-expanded target lists.

`python verifier/mass_denominators.py` evaluates the recorded cofactor and split-tree denominators at (3,7). The scan finds 506 distinct vanishing denominator expressions across ten files. These identities cannot specialize there as written; the scan does not certify a named instance of the theorem.

The trusted arithmetic includes Python, the verifier code, FLINT/python-flint and SymPy for the companion identities. The completeness audit additionally uses cddlib/GMP, Z3, and Yices. Singular and gfan generate evidence but are not called during replay. These checks are not a proof-assistant formalization.

The second solver replay uses standard SMT-LIB QF_LRA and a Boolean cardinality encoding, independent of Z3's pseudo-Boolean extension:

```sh
python verifier/fan_second_solver.py
```

It checks all six normalized symmetry cases, and validates complete and missing-ray tropical-line controls. `UNKNOWN`, timeouts and partial runs fail. The input half-spaces are regenerated in a full replay.

All six cases returned `UNSAT` in both Z3 5.1.0 and Yices 2.7.0. The second-solver results and query hashes are in `data/fan_second_solver_audit.json`, with its standard SMT-LIB base query in `data/fan_second_solver.smt2.gz`.

## Explicit mass-factor inventory

`python verifier/mass_factor_inventory.py` factors all archived cofactor/split denominators and the fourteen companion batches at `c=1`, checks factorization reconstruction, and records exact positive-zero witnesses or opposite-sign endpoints. Its output defaults to `logs/`; the checked inventory is `data/mass_factor_inventory.json`.

The cofactor denominators yield 798 irreducible factors: 295 have coefficients of one sign and cannot vanish at positive masses, 421 have certified positive zeros, and 82 remain unclassified. These are conditions for the recorded identities, not a classification of masses with infinitely many configurations. Alternative identities can remove some conditions. The modular and elimination spreading conditions have not all been extracted, so this is not the complete exceptional polynomial.

The rational pairs `(7,13)` and `(11,17)` avoid all scanned cofactor and companion factors. They remain candidates: named instances additionally require point-specific structural, elimination, and modular checks.

The stable filenames `verify_lemS.py` and `verify_lemS_rep.py` refer to the modular-star criterion; their names retain the historical lemma label.
