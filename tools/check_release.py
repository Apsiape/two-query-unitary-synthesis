"""Bounded, read-only repository consistency checks, not theorem verification."""
from pathlib import Path
import re
import subprocess
import sys
import json

root = Path(__file__).resolve().parents[1]
paths = ("README.md", "paper/paper.md", "paper/paper.tex",
         "notes/garbage-to-clean.md", "notes/permutation-transport.md",
         "STATUS.md", "CITATION.cff", "lean/README.md")
text = {}
for name in paths:
    data = (root / name).read_bytes()
    assert len(data) < 1024 * 1024, name
    s = data.decode("utf-8")
    assert all(ord(c) >= 32 or c in "\r\n\t" for c in s), (name, "control byte")
    text[name] = s
version = re.search(r"^version: ([0-9.]+)$", text["CITATION.cff"], re.M).group(1)
for name in ("README.md", "paper/paper.md", "paper/paper.tex"):
    assert version in text[name], (name, "version mismatch")
metadata = json.loads((root / "publication/zenodo-metadata.json").read_text(encoding="utf-8"))
assert metadata["version"] == version
assert metadata["upload_type"] == "publication" and metadata["publication_type"] == "preprint"
assert metadata["creators"] == [{"name": "Douglas, Seth"}]
assert metadata["access_right"] == "open" and metadata["license"] == "cc-by-4.0"
assert "doi" not in metadata, "Do not invent or reuse an unverified DOI"
md = text["paper/paper.md"]
for n in range(1, 10):
    assert len(re.findall(r"^## " + str(n) + r"\. ", md, re.M)) == 1, (n, "heading count")
tex = text["paper/paper.tex"]
assert "\\label{app:transport}" in tex and "\\label{app:pair}" in tex
assert not re.search(r"\\(?:input|include|bibliography)\s*\{", tex)
labels = re.findall(r"\\label\{([^}]+)\}", tex)
assert len(labels) == len(set(labels)), "duplicate label"
assert set(re.findall(r"\\(?:eqref|ref)\{([^}]+)\}", tex)) <= set(labels), "missing label"
assert "Every proof was re-derived by the author and is machine-checked" not in tex
assert "no width argument that is linear" not in tex
assert "C_S\\le2.09" not in tex
assert "Packing is over $g$ alone" in tex, "lost fixed-architecture quantifier"
for name in ("paper/paper.md", "notes/garbage-to-clean.md"):
    assert "Nothing below degree 4 can suffice" not in text[name]
mdfiles = [root / p for p in paths if p.endswith(".md")]
subprocess.run([sys.executable, str(root / "verify/md_math_lint.py"),
                *map(str, mdfiles)], check=True, timeout=15)
print("PASS: version, unique headings/labels, self-contained TeX, scope and Markdown lint")
