"""Install the published certificate archive after checking all file hashes."""
from pathlib import Path
import argparse, hashlib, json, os, shutil, tarfile, tempfile, urllib.request

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = "brandonlign/tropcc"
RELEASE = "v1.3.4"

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Use a local archive instead of downloading")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "certificate_manifest.json").read_text())
    with tempfile.TemporaryDirectory(prefix="tropcc-certificates-") as tmp:
        archive = args.archive
        if archive is None:
            archive = Path(tmp) / manifest["archive"]
            print("Downloading certificate archive...", flush=True)
            url = ("https://github.com/" + REPOSITORY + "/releases/download/"
                   + RELEASE + "/" + manifest["archive"])
            urllib.request.urlretrieve(url, archive)
        if digest(archive) != manifest["sha256"]:
            raise ValueError("Archive SHA-256 mismatch")
        expected = manifest["files"]
        found = set()
        with tarfile.open(archive, "r:xz") as tar:
            for member in tar:
                if not member.isfile() or member.name not in expected or member.name in found:
                    raise ValueError("Unexpected archive member: " + member.name)
                target = (ROOT / member.name).resolve()
                if not target.is_relative_to(ROOT):
                    raise ValueError("Unsafe archive path")
                target.parent.mkdir(parents=True, exist_ok=True)
                stream = tar.extractfile(member)
                with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as out:
                    temporary = Path(out.name)
                    shutil.copyfileobj(stream, out)
                if digest(temporary) != expected[member.name]:
                    temporary.unlink()
                    raise ValueError("Certificate SHA-256 mismatch: " + member.name)
                os.replace(temporary, target)
                found.add(member.name)
        if found != set(expected):
            raise ValueError("Archive is incomplete")
    print("Installed and verified", len(found), "certificate files.")

if __name__ == "__main__":
    main()
