# Finiteness of five-body central configurations for generic equal-pair masses

**Brandon Li · Ethan Yeroushalmi**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23197573.svg)](https://doi.org/10.5281/zenodo.23197573)

**Paper:** [`paper/tropcc.pdf`](paper/tropcc.pdf) · [LaTeX source](paper/tropcc.tex) · [Certificate map](CERTIFICATES.md)

## The result

Smale's sixth problem asks whether every choice of positive masses gives only finitely many central configurations. Albouy and Kaloshin (Annals, 2012) answered this for five bodies in the plane, except for an exceptional set of masses. One part of that set is the equal-pair family, with masses `(a, a, b, b, c)`.

**Theorem.** Outside a proper algebraic hypersurface of positive triples `(a, b, c)`, the planar five-body problem with masses `(a, a, b, b, c)` has finitely many central configurations up to similarity.

The theorem does not cover every mass in the family. The exceptional hypersurface is proved to exist but is not written down explicitly, and the configurations are not counted.

## How the proof works

1. **Tropical fan.** If a mass choice had infinitely many configurations, the growth rates of the ten mutual distances along some curve of solutions would give a weight vector. That vector must lie in a polyhedral fan of 16,438 cones.
2. **Exclude cones.** For 16,340 cones, an exact certificate shows that the leading-order equations have no solution with all distances nonzero. The certificates are Laurent unit identities, eliminations, split trees, and a modular star criterion.
3. **Half-space.** The remaining 98 cones (77 rays) all have positive pairing with `η = (-23, 1, 1, 51, 1, 1, 51, 24, 13, 13)`, so the solution variety is finite.
4. **Physics.** This is done at a fixed potential. Real configurations have only finitely many potential values, which gives finitely many configurations overall.

## Verifying the proof

Requirements:

- Python 3.10 or later
- cddlib with GMP (`cddexec_gmp`)
- Yices 2.7 (`yices-smt2`)
- On macOS: `brew install cddlib yices2`
- On Debian/Ubuntu: `libcdd-tools python3-z3`, plus Yices from its official distribution

Singular and gfan are used to *generate* the evidence but are not needed to check it.

```sh
python3 -m venv --system-site-packages .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/fetch_certificates.py   # downloads and hash-checks the certificate archive (~1.3 GB extracted)
python scripts/verify_all.py           # full replay
```

The certificate archive (167 MB) is attached to the [v1.3.9 release](https://github.com/brandonlign/tropcc/releases/tag/v1.3.9). `fetch_certificates.py` downloads it, checks its SHA-256 digest and checks every extracted file against `certificate_manifest.json`. For an offline copy, pass `--archive /path/to/certificates.tar.xz`.

A successful run rechecks every certificate and every symmetry transport. It then rebuilds the coverage from scratch and confirms the 98-cone residual and the half-space pairing. Use the full replay, not `coverage.py --fast`, which trusts cached markers.

### Tested environment

The full replay was rerun on October 8, 2026, at v1.3.8 (commit `eb2645f`) and passed all 15 stages with exit status 0 and `ALL PASS`. See the [full replay summary](verification_logs/full_replay_v1.3.8.log) and [detailed stage logs](verification_logs/full_replay_v1.3.8_details.log). All 114 source hashes and 1,865 certificate hashes matched before and after the run. This replay used the existing machine; the earlier fresh-clone run is recorded in [the v1.3.1 log](verification_logs/full_replay_v1.3.1.log).

| | |
|---|---|
| Machine | Apple M3, 16 GB RAM, macOS 26.6 |
| Software | Python 3.14.6, python-flint 0.9.0, SymPy 1.14.0, Z3 5.1.0, Yices 2.7.0, cddlib 0.94n |
| Runtime | 28 minutes 17 seconds, single-threaded |
| Peak memory | about 3.8 GB |
| Disk | about 1.3 GB after extracting the certificates |

The [current manuscript (v1.3.10)](https://doi.org/10.5281/zenodo.23249007) cites the verified [v1.3.9 package](https://doi.org/10.5281/zenodo.23248651). Its equations, certificates, fan, and computational checkers are unchanged from the recorded full v1.3.8 replay.

## Repository layout

| Path | Contents |
|---|---|
| `paper/` | Manuscript (PDF and LaTeX) |
| `data/` | Equations, support fan, audits, and separating vectors |
| `verifier/` | Exact checks for certificates, fan, symmetries, and coverage |
| `companion/` | The 10,965-cone Laurent ledger and its checkers |
| `scripts/` | Certificate download and full replay |
| `verification_logs/` | Log of the recorded full replay |
| `CERTIFICATES.md` | Map from each claim in the paper to its evidence and checker |

## What is trusted

The replay uses exact arithmetic throughout: Python, FLINT (python-flint), and SymPy. It is not a proof-assistant formalization. The fan is checked in two ways:

- `verifier/fan_structure.py` rebuilds all 16,438 cones and 44,998 facet incidences without any external solver.
- The support-completeness audit uses cddlib/GMP together with two independent SMT solvers, Z3 5.1.0 and Yices 2.7.0. All six symmetry cases returned `UNSAT` in both solvers. A timeout or `UNKNOWN` counts as a failure.

## Known limits

- **Not every mass is covered.** Some exclusions come from elimination or the modular criterion. For those, the mass polynomials to avoid are proved to exist but are not extracted.
- **No specific mass is certified.** `verifier/mass_factor_inventory.py` lists the explicit denominator factors, and the pairs `(7,13)` and `(11,17)` avoid all of them. They are still only candidates, not certified instances of the theorem.
- **Conditional transfers have their own exceptions.** The exceptional hypersurface includes the cone-325 obstruction `(16a²-b²)(4a²+b²)` and its image under `a ↔ b`. For positive masses these require avoiding `b=4a` and `a=4b`. The all-positive exclusion for cone 346 does not remove these transfer conditions. The factor inventory records both full complex obstructions; the full replay regenerates the inventory under `logs/mass_factor_inventory.json`.
- **An old file label.** The scripts `verify_lemS*.py` implement the modular star criterion; the names keep an old lemma label.

## Citation

The paper cites the verified [v1.3.9 package](https://doi.org/10.5281/zenodo.23248651). The [concept record](https://doi.org/10.5281/zenodo.23197573) covers all versions; `CITATION.cff` intentionally cites that record without a version number. Code is released under the MIT License; the paper under CC BY 4.0.

Contact: brandon.li.gn@gmail.com
