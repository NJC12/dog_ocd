# s_het comparison of ST6 dog and human genes

All 87 Dog and 48 Human rows in `NATURE_SupplementaryTables.xlsx`, sheet ST6,
were included. The labels identify the discovery gene set; s_het is a human
estimate in both groups. The 78-gene orthology restriction used by the earlier
GWAS enrichment analysis was not applied. The four ST6 non-one-to-one genes
and five previously unconfirmed dog orthologs were retained by their named
human counterparts. This compares the full supplementary gene lists, rather
than establishing a new one-to-one orthology set.

## Posterior means and mapping

Inputs are the six `posterior_densities.*.txt.gz` files and
`elife-83172-supp2-v2.txt` in `~/data/constraint/shet/agarwal/` from
[Agarwal et al.](https://elifesciences.org/articles/83172).
The summary contains MAP estimates; posterior means were calculated from the
densities. The [authors' calculation code](https://github.com/agarwal-i/loss-of-function-fitness-effects/blob/ae3ffcbd947c578c5b23eb3db81076877ecbad5a/code/calculate_posterior_densities.R)
defines these densities on x = log10(s_het).

The approximate posterior mean is ∫10^x f(x) dx / ∫f(x) dx, calculated by
trapezoidal integration on the supplied 1,024-point grid from -7 to 0.
The y-axis is log10 of this arithmetic posterior mean. The density integral
was checked and normalized separately for each gene; observed integrals range
from 0.994132 to 1.000591. Estimates reflect the published smoothed, truncated
density grids, rather than original posterior draws.

Five renamed symbols were matched through the exact Ensembl transcript ID in
the Agarwal summary and GENCODE v50: ATP5MC1/ATP5G1 (ENST00000393366),
H4C11/HIST1H4J (ENST00000355057), KATNIP/KIAA0556 (ENST00000261588),
MAIP1/C2orf47 (ENST00000295079), and TENT4A/PAPD7 (ENST00000230859).
[ADGRG1/GPR56](https://www.ncbi.nlm.nih.gov/gene/9289) and
[PROSER3/C19orf55](https://www.ncbi.nlm.nih.gov/gene/148137) were checked against
NCBI Gene. HLA-B is stored as `HLA.B` in the posterior file. FAM155A maps to
NALF1 in the GWAS table using the existing orthology audit; its Agarwal name
remains FAM155A. COL6A5 is absent from both the summary and posterior files;
it is retained with missing s_het in the output and excluded from the test
and scatter plot. All other 134 genes have finite, positive autosomal estimates.

## Ranking and comparison

All 134 available genes were ranked together by posterior mean s_het,
descending (rank 1 highest; average ranks for ties). The one requested
comparison uses a two-sided, asymptotic Mann–Whitney U test with tie and
continuity corrections, with dog genes as the first sample. No multiple-test
correction is applied to this single test. A separate rank-sum identity
confirmed the U statistic.

| ST6 group | Genes tested | Median posterior mean s_het | Median pooled rank |
| --- | ---: | ---: | ---: |
| Dog | 86 | 0.0415423 | 64.5 |
| Human | 48 | 0.0114333 | 78.5 |

U_dog = 2,268; two-sided P = 0.345007. U/(86×48) = 0.549419 estimates the
probability that a randomly selected dog-list gene has higher s_het than a
randomly selected human-list gene, counting ties as one half. This test does
not show a significant distributional difference. It compares point estimates
and does not propagate their posterior uncertainty.

## Figure and outputs

The scatter plot contains all 134 tested genes, colored by ST6 group. Blue
and orange squares mark the dog and human group medians, respectively.
Each square uses the median raw window p-value and median posterior mean
s_het among that group's plotted genes, transformed to the displayed axes.
Its x-axis is -log10 of `min_p_500kb` in the existing `per_gene_ocd.csv`: the
minimum OCD SNP p-value within the inclusive ±500 kb interval around each
GENCODE gene-level TSS. All 135 original ST6 genes have these window values.
Overlapping windows can share p-values, and the human list was selected using
human OCD results; the scatter plot is descriptive.

Dog points use dark blue for ST6 `dogCD significance = significant` and light
blue for `suggestive`; human points remain orange. These categories come from
the dog GWAS, not the human OCD window p-value on the x-axis. The dog median
square continues to summarize both dog categories together.

The pooled best-fit line is ordinary least squares with an intercept, fitting
log10(posterior mean s_het) against -log10(minimum OCD window p-value) across
all 134 plotted genes, with equal weight per gene. Median squares are summaries
and are not included as additional observations. Its slope is printed on the
plot; full-precision slope, intercept, R², and sample size are saved in
`shet_pooled_regression.csv`. This line is descriptive and does not account for
shared windows, posterior uncertainty, or gene-list selection.

The upper marginal panel shows density against -log10(window p-value), sharing
the scatter plot's x-axis. The right marginal panel shows density against
log10(posterior mean s_het), sharing its y-axis. Each compares the same 86 dog
genes (significant and suggestive combined, blue) with the same 48 human genes
(orange). Densities are Gaussian kernel estimates fitted separately to each
group in these transformed coordinates, using SciPy's Scott bandwidth rule
and equal weight per gene. Each KDE has unit total mass before display limits;
these are distributions of point estimates, not pooled posterior densities.
Axes extend to accommodate three kernel bandwidths beyond the observations,
limited to -log10(P) ≥ 0 and log10(s_het) ≤ 0. Evaluated density curves are
saved in `shet_marginal_densities.csv`.

Run `.venv/bin/python src/compare_shet.py` from the project root. Section 5 of
`dog_ocd_analysis.ipynb` runs the same analysis and displays its results.

1. `shet_gene_ranks.csv`: all 135 genes, mapping evidence, estimates, ranks,
   window p-values, and missingness status.
2. `shet_mann_whitney.csv`: full-precision test results and group summaries.
3. `shet_vs_ocd_500kb.png` and `.svg`: scatter plot.
4. `shet_audit.json`: input SHA-256 hashes, density checks, results, and
   Python/package versions; `analysis_output.txt`: printed full ranking.

The analytical integration test checks a uniform log-density against its
closed-form linear-scale mean. Further tests cover invalid densities, missing
estimates, ties, and U orientation. The full 10-test project suite passed.
The rendered plot was visually inspected and opened in the default image app.
