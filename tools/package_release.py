"""Package the inspected working tree, without committing, pushing or publishing.

The source ZIP contains the explicit public source roots only. Its manifest says
whether the files differ from HEAD; a working-tree archive is never called a commit.
Run check_release.py and build_paper.py before this script.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import zipfile

root = Path(__file__).resolve().parents[1]
out = root / "release"
out.mkdir(exist_ok=True)
version = re.search(r"^version: ([0-9.]+)$", (root / "CITATION.cff").read_text(), re.M).group(1)
top_files = {".gitattributes", ".gitignore", ".zenodo.json", "README.md", "CHANGELOG.md", "CITATION.cff",
             "LICENSE", "LICENSE-CODE", "STATUS.md", "VALIDATION.md", "requirements.txt"}
source_roots = {"paper", "notes", "verify", "lean", "tools", "audit", "publication"}
listed = subprocess.check_output(
    ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
    cwd=root, timeout=15).decode().split("\0")
paths = []
for name in sorted(set(filter(None, listed))):
    p = Path(name)
    if name in top_files or p.parts[0] in source_roots:
        assert not any(x in p.parts for x in (".lake", "__pycache__", "..")), name
        path = root / p
        assert path.is_file() and not path.is_symlink(), name
        assert path.stat().st_size < 10 * 1024 * 1024, name
        paths.append(name)
assert "audit/2026-09-13-independent-audit.md" in paths, "Audit report missing"
assert "publication/RELEASE-CHECKLIST.md" in paths, "Release checklist missing"
archive = out / f"two-query-unitary-synthesis-v{version}.zip"
manifest = {"version": version, "kind": "working-tree snapshot",
            "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root,
                                                     timeout=15).decode().strip(),
            "working_tree_dirty": bool(subprocess.check_output(
                ["git", "status", "--porcelain"], cwd=root, timeout=15).strip()), "files": {}}
with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
    for name in paths:
        data = (root / name).read_bytes()
        manifest["files"][name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
        info = zipfile.ZipInfo(f"two-query-unitary-synthesis-v{version}/{name}", (2026, 9, 13, 0, 0, 0))
        info.compress_type = zipfile.ZIP_DEFLATED
        z.writestr(info, data)
manifest_path = out / "SOURCE-MANIFEST.json"
manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
artifacts = [root / "paper/paper.pdf", out / "arxiv-source.zip", archive, manifest_path]
assert all(p.is_file() for p in artifacts)
checksums = "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(root).as_posix()}\n"
                    for p in artifacts)
(out / "SHA256SUMS.txt").write_text(checksums, encoding="utf-8")
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
print(f"PASS: {len(paths)} public source files; ZIP integrity; SHA-256 manifest")
print(archive)
