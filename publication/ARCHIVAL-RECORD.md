# Version 1.0.4 archival publication

Verified on 2026-09-13 (America/Los_Angeles; release timestamp
2026-09-14T04:09:18Z).

- Version DOI: [10.5281/zenodo.22740995](https://doi.org/10.5281/zenodo.22740995).
- All-versions DOI: [10.5281/zenodo.22740994](https://doi.org/10.5281/zenodo.22740994).
- GitHub release: [v1.0.4](https://github.com/Apsiape/two-query-unitary-synthesis/releases/tag/v1.0.4).
- Released commit: `99195f88f4b562958b629b48334e5fda71108c6e`.
- Zenodo record type: preprint; author: Seth Douglas; version: 1.0.4.

## Checks completed

The public Zenodo API returned the correct title, author, version, preprint type,
and explicit file-specific license description. The downloaded archive's MD5
matched Zenodo's advertised checksum. All 46 archived files matched the released
Git commit byte for byte, including the audited PDF and TeX.

The five GitHub release assets (PDF, standalone LaTeX ZIP, full source ZIP,
source manifest, and checksums) matched their local SHA-256 digests. The local
source ZIP is a separately packaged Windows checkout: seven text files differ
from Git's canonical blobs only by CRLF versus LF. Its manifest describes that
local package, not Zenodo's automatically generated ZIP. No substantive contents
differ.

Audited manuscript hashes:

```text
paper/paper.tex  b9bf7d21198a032e270a1322f4a590268c523975d03cc7ceffffe44f7a819bb8
paper/paper.pdf  1529c00e91c5153b8f03fbaef6283f2c149ce9294574170faa5f45bed79e508d
```

The repository-level `.zenodo.json` overrides the integration's limited CFF
importer. The record-level CC BY 4.0 label describes the manuscript/prose;
software remains MIT, explicitly stated in both the description and notes.
No file was relicensed.

The Zenodo integration archived the repository ZIP, which includes the PDF; it
did not separately ingest the GitHub release attachments. The PDF is also
directly downloadable from the GitHub release.

## Post-release citation update

README, CITATION.cff, and publication records on `main` now link the verified DOI.
These citation-only updates do not move the release tag, replace release assets,
change the manuscript, or create another Zenodo version. The archived snapshot
necessarily precedes the DOI-link update.

This is public archival publication, not external mathematical peer review,
journal acceptance, or a resolution of full AK. See VALIDATION.md for proof and
verification scope. No mathematical test suite was rerun merely for DOI links.
