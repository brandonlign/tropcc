# Certificate map

| Claim | Evidence | Check |
|---|---|---|
| Polynomial system and monomial supports | `data/sys.json`, `data/supp.in` | `suppcheck.py`, companion `verify_equations.py` |
| Independent support-fan agreement | `data/fan.out`, `verifier/data/fan062.out.gz` | `fancheck.py` |
| Pair swaps and mass-exchange symmetry | `verifier/orbits.json`, `verifier/sigma.json` | `symcheck.py`, `symcheck2.py` |
| 10,965 companion exclusions | `companion/results/` | `verify_companion.py` |
| Laurent and resultant identities | `data/exact/{cert,dcert,gcert,ecert}/` | `verify.py`, `verify_res.py` |
| Split-tree identities | `data/exact/split/` | `verify_split.py` |
| Elimination and symmetry transfers | Initial forms reconstructed from the equations and fan | `verify_334.py`, `verify_49.py`, `elim_sigma.py` |
| Modular specialization for cells 81, 121, 206 | `data/lemS/`, including exceptional-branch identities | `verify_lemS.py`, `verify_lemS_rep.py`, `b0_verify.py`, `cluster_mod.py` |
| Complete coverage and residual half-space | `data/balance/`, `verifier/coverage_report.json` | `coverage.py` |

Verifier paths in the last column are relative to `verifier/` unless stated otherwise. Large evidence files are installed from the release archive by `scripts/fetch_certificates.py`.

The complete replay is `python scripts/verify_all.py`. It rechecks the companion proofs, every modular star, characteristic-zero identities, elimination branches, and residual ray pairings. Run the full check; `coverage.py --fast` trusts convenience markers and is not a proof check. Replay logs are written locally under `logs/`.
