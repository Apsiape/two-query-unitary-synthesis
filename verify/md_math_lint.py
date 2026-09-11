r"""md_math_lint.py - GitHub math lint. Rules learned in blood:
GitHub markdown processes emphasis (*, _) INSIDE $$ blocks and $...$
spans, eating asterisks/underscores before the math renderer runs.
Therefore: display math MUST use ```math fences; inline math MUST use
$`...`$ backtick-protected spans. Also: \{ \} lose their backslashes
(use \lbrace/\rbrace), and \! adjacent to delimiters errors."""
import re, sys
bad = 0
for path in sys.argv[1:]:
    text = open(path, encoding="utf-8").read()
    if "$$" in text:
        print(f"{path}: bare $$ display math - use ```math fences"); bad += 1
    for i, l in enumerate(text.splitlines(), 1):
        if re.search(r"(?<![`$])\$(?![`$])[^$]*(?<![`$])\$(?![`$])", l) and "```" not in l:
            print(f"{path}:{i}: unprotected inline $...$ - use $`...`$"); bad += 1
        if r"\!\left" in l or r"\!\right" in l or r"\!\Big" in l:
            print(f"{path}:{i}: \\! adjacent to delimiter"); bad += 1
        if re.search(r"\\{|\\}", l) and "lbrace" not in l and "rbrace" not in l and "```" not in l:
            pass  # inside ```math fences \{ \} are safe; only flag in inline spans
    for block in re.findall(r"```math(.*?)```", text, re.S):
        if block.count(r"\left") != block.count(r"\right"):
            print(f"{path}: unbalanced left/right in a math fence"); bad += 1
sys.exit(1 if bad else 0)
