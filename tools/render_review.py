"""Render all PDF pages and make labelled contact sheets for human visual review."""
from pathlib import Path
import subprocess
import re
from PIL import Image, ImageDraw

root = Path(__file__).resolve().parents[1]
version = re.search(r"^version: ([0-9.]+)$", (root / "CITATION.cff").read_text(), re.M).group(1)
out = root / "tmp" / "pdfs" / version
out.mkdir(parents=True, exist_ok=True)
subprocess.run(["pdftoppm", "-r", "85", "-png", str(root / "paper/paper.pdf"),
                str(out / "page")], check=True, timeout=60, capture_output=True)
pages = sorted(out.glob("page-*.png"))
assert pages
for start in range(0, len(pages), 6):
    batch = pages[start:start+6]
    canvas = Image.new("RGB", (1800, 1760), "gray")
    draw = ImageDraw.Draw(canvas)
    for n, path in enumerate(batch):
        with Image.open(path) as page:
            page.thumbnail((590, 830))
            x, y = (n % 3) * 600, (n // 3) * 880
            canvas.paste(page, (x+5, y+25))
            draw.text((x+8, y+8), path.stem, fill="white")
    dest = out / f"contact-{start//6+1}.png"
    canvas.save(dest)
    print(dest)
print(f"Rendered {len(pages)} pages; contact sheets do not replace visual inspection.")
