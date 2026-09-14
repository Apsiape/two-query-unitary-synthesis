"""Build the self-contained arXiv source in an isolated temporary directory.

Uses installed pdflatex; does not install packages, release, or upload anything.
Invoke under the campaign memory/time wrapper on the local Windows workstation.
"""
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / "paper" / "paper.tex"
dest = ROOT / "release"
dest.mkdir(exist_ok=True)
with tempfile.TemporaryDirectory(prefix="q2-paper-") as tmp:
    work = Path(tmp)
    shutil.copy2(source, work / "paper.tex")
    for _ in range(3):
        run = subprocess.run(
            ["pdflatex", "-no-shell-escape", "-interaction=nonstopmode",
             "-halt-on-error", "paper.tex"],
            cwd=work, capture_output=True, text=True, errors="replace", timeout=40)
        if run.returncode:
            raise RuntimeError(run.stdout[-6000:] + run.stderr[-1000:])
    log = (work / "paper.log").read_text(errors="replace")
    forbidden = ("Undefined control sequence", "undefined references",
                 "undefined on input line", "Overfull \\hbox", "Overfull \\vbox",
                 "multiply defined")
    assert not any(s in log for s in forbidden), "\n".join(
        line for line in log.splitlines() if any(s in line for s in forbidden))
    shutil.copy2(work / "paper.pdf", ROOT / "paper" / "paper.pdf")
    shutil.copy2(work / "paper.log", dest / "build.log")
with zipfile.ZipFile(dest / "arxiv-source.zip", "w", zipfile.ZIP_DEFLATED) as z:
    z.write(source, "paper.tex")
print("PASS: isolated three-pass PDF build; no missing references or overfull boxes")
print("Source bundle: release/arxiv-source.zip (paper.tex only; bibliography embedded)")
