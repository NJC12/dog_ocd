# Input publication notes

## Supplementary workbook

The gene-set source is Tengvall et al., *Dog-human gene network colocalization
links synaptic signalling and ubiquitination to compulsive disorders*, bioRxiv
version 1, 2026-10-12, DOI 10.64898/2026.10.12.758413.

[Public supplementary material](https://www.biorxiv.org/content/10.64898/2026.10.12.758413v1.supplementary-material)

The project workbook was compared with the publicly downloaded `media-2(2).xlsx`
on 2026-11-16. Both have the same index and 26 supplementary-table sheets.
ST5 and ST6, used by this analysis, have identical cell contents. Across the
workbook, 12 cells differ: seven heading edits, four equivalent constant numeric
formulas, and one floating-point representation difference. File checksums
differ; no byte-identity claim is made. The original project workbook remains
unchanged. The comparison report is in
`data/generated/supplementary_public_comparison.json`.

The accompanying `media-1-3.pdf` contains supplementary methods for dog brain
regulatory-element profiling, compulsive-behaviour factors and related study
methods. It was inspected as source context and is not an analysis input.

The bioRxiv metadata record identifies the preprint license as CC BY-NC-ND 4.0.
This repository does not relicense the supplementary material. Preserve creator
attribution and consult the source's license when reusing it.

## All by All burden data: excluded

Associations can be viewed at https://allbyall.researchallofus.org/ by searching
for phenotype code `obsessive_compulsive_disorder`, selecting ancestry `meta`,
predicted loss-of-function (pLoF), and MAF < 0.001.

The [official FAQ](https://allbyall.researchallofus.org/faq), also available in the
[browser's source repository](https://github.com/broadinstitute/all-by-all-aou-browser/blob/main/frontend/src/client/about/faq.md),
states that downloads of public browser results are not permitted. It describes
a separate route for registered researchers to download publication summary
tables through the Researcher Workbench. The same wording was present in the
deployed browser JavaScript when checked on 2026-11-16.

The browser-exported CSV and all generated `notebook/results/burden/` outputs
are excluded from the public repository. This includes copied source association
statistics in `gene_ranks.csv`, plots and aggregate analysis results. The local
files remain intact. The analysis code and source-provenance TXT are included;
the main notebook does not use these burden data.

The public branch starts with a new root commit so the excluded data are absent
from its entire history. Earlier local commits are retained on the local `main`
branch and must not be pushed because they contain the excluded files. Push only
the `public-main` branch to GitHub's `main` branch.
