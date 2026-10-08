"""Replay the complete proof from the distributed equations and certificates."""
from pathlib import Path
import gzip, json, subprocess, sys

ROOT = Path(__file__).resolve().parents[1]

def main():
    manifest = json.loads((ROOT / "certificate_manifest.json").read_text())
    missing = [name for name in manifest["files"] if not (ROOT / name).is_file()]
    if missing:
        raise SystemExit("Install certificates first: python scripts/fetch_certificates.py")
    logs = ROOT / "logs"
    logs.mkdir(exist_ok=True)
    independent = ROOT / "verifier/data/fan062.out"
    with gzip.open(str(independent) + ".gz", "rb") as stream:
        independent.write_bytes(stream.read())
    commands = [
        ("review_controls", ["verifier/review_checks.py"]),
        ("supports", ["verifier/suppcheck.py"]),
        ("pair_symmetry", ["verifier/symcheck.py"]),
        ("mass_symmetry", ["verifier/symcheck2.py"]),
        ("fan", ["verifier/fancheck.py", str(independent)]),
        ("fan_soundness", ["verifier/fan_soundness.py"]),
        ("fan_completeness", ["verifier/fan_completeness.py", "--symmetry-break"]),
        ("companion", ["verifier/verify_companion.py"]),
        ("modular_stars", ["verifier/verify_lemS.py", "81", "121", "206"]),
        ("coverage", ["verifier/coverage.py"]),
        ("cone_accounting", ["verifier/cone_accounting.py"]),
    ]
    for name, command in commands:
        print("Checking", name, "...", flush=True)
        with (logs / (name + ".log")).open("w") as out:
            result = subprocess.run([sys.executable] + command, cwd=ROOT, stdout=out, stderr=subprocess.STDOUT)
        if result.returncode:
            raise SystemExit("FAILED " + name + "; see logs/" + name + ".log")
        print("PASS", name, flush=True)
    report = json.loads((ROOT / "verifier/coverage_report.json").read_text())
    assert report["mode"] == "full" and report["pointed"]
    assert report["total_cells"] == 16438 and len(report["residual_cells"]) == 98
    assert len(report["residual_rays"]) == 77
    assert report["c"] == [-23,1,1,51,1,1,51,24,13,13]
    print("ALL PASS: 16,340 excluded cones; 98 residual cones; 77 rays; exact positive pairings.")

if __name__ == "__main__":
    main()
