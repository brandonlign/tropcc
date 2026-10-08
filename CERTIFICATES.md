# Certificate map

| Claim | Evidence | Check |
|---|---|---|
| Polynomial system and monomial supports | `data/sys.json`, `data/supp.in` | `suppcheck.py`, companion `verify_equations.py` |
| Cross-version support-fan agreement | `data/fan.out`, `verifier/data/fan062.out.gz` | `fancheck.py` |
| Independent support-fan containment in both directions | `data/fan_completeness_audit.json`, checked rational query | `fan_soundness.py`, `fan_completeness.py` |
| Second solver completeness replay | `data/fan_second_solver_audit.json`, standard SMT-LIB base query | `fan_second_solver.py` (Yices) |
| Exact Newton normal regions, facets, fan property, and full fan symmetries | `data/fan_structure_audit.json` | `fan_structure.py` |
| Recorded mass denominators at (3,7) | `data/mass_denominator_audit.json` | `mass_denominators.py` |
| Disjoint proof-type accounting | `data/cone_accounting.json` | `cone_accounting.py` |
| Pair swaps and mass-exchange symmetry | `verifier/orbits.json`, `verifier/sigma.json` | `symcheck.py`, `symcheck2.py` |
| 10,965 companion exclusions | `companion/results/` | `verify_companion.py` |
| Laurent and resultant identities | `data/exact/{cert,dcert,gcert,ecert}/` | `verify.py`, `verify_res.py` |
| Split-tree identities | `data/exact/split/` | `verify_split.py` |
| Elimination and symmetry transfers | Initial forms reconstructed from the equations and fan | `verify_334.py`, `verify_49.py`, `elim_sigma.py` |
| Elimination inputs, mass reduction, and specialization counterexamples | `data/elimination_forms.json` | `review_checks.py` |
| Modular specialization for cells 81, 121, 206 | `data/lemS/`, including exceptional-branch identities | `verify_lemS.py`, `verify_lemS_rep.py`, `b0_verify.py`, `cluster_mod.py` |
| Complete coverage and residual half-space | `data/balance/`, `verifier/coverage_report.json` | `coverage.py` |

Verifier paths in the last column are relative to `verifier/` unless stated otherwise. Large evidence files are installed from the release archive by `scripts/fetch_certificates.py`.

The complete replay is `python scripts/verify_all.py`. It rechecks the companion proofs, every modular star, characteristic-zero identities, elimination branches, and residual ray pairings. Run the full check; `coverage.py --fast` trusts convenience markers and is not a proof check. Replay logs are written locally under `logs/`.
