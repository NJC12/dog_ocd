# GWAS enrichment using ST6 TOGA orthologs

## Gene selection and scope

The analysis selects Dog rows from `NATURE_SupplementaryTables.xlsx`, sheet ST6,
whose `TOGA (Human:Dog)**` value is `1:1` after trimming whitespace. ST6 contains
87 unique Dog genes: 83 labeled `1:1` and four labeled `Not 1:1`. The selected
83 comprise 43 dog-GWAS-significant and 40 dog-GWAS-suggestive genes.

Human scoring windows are identified by matching the selected symbols to the
autosomal protein-coding GENCODE v50 annotation. The former symbol `FAM155A`
is mapped to `NALF1` (`ENSG00000204442`). Ensembl homology classifications do
not filter the TOGA selection. This restores `EMC2`, `PPA2`, `SH3GL3`,
`PLXNA2` and `PROSER3` relative to the earlier 78-gene analysis. `EMC2`, `PPA2`,
`PLXNA2` and `PROSER3` are dog-GWAS-significant; `SH3GL3` is suggestive. The selected symbols, mapped
human IDs and significance labels are exported in `toga_dog_gene_set.csv`.

The analysis covers human OCD, CAD, BIP, ADHD, MDD and SCZ GWAS at inclusive
±1 Mb and ±500 kb windows around the gene-level TSS. Established genome-wide
window scores are reused after checking their provenance and annotation;
dog membership labels are replaced before tests and plotting. Window p-values,
SNP counts, lead SNPs and gene coordinates retain their established definitions.
The source data and coordinate processing are documented in
[the original methods](../../METHODS.md) and the
[project README](../../../../README.md).

## Inference and figures

For each trait and window, covered genes are ranked by their minimum window
p-value, with rank 1 strongest and average ranks for ties. Three Dog sets are
tested: significant, suggestive and combined. All use the same reference group
of other covered genes, excluding the entire 83-gene Dog set. One-sided
asymptotic Mann–Whitney tests ask whether Dog genes have stronger GWAS signal.
Holm adjustment covers only the six OCD tests: two windows × three Dog sets.
CAD, BIP, ADHD, MDD and SCZ are control traits and report raw p-values without
Holm adjustment. Their Holm values are blank in CSV exports and shown as a
dash in notebook/HTML tables. Repeated OCD displays are the same six tests.

The Human-set analysis retains all 48 ST6 Human genes and compares them with
all other covered genes. Its Holm family is separate and contains 12 tests:
six traits × two windows. Changes in Dog membership do not alter the Human-set
definition or its background.

Circular-shift tests use five SNP-count strata within each chromosome and
percentile ranks within chromosome/stratum. Labels are rotated over TSS-ordered
covered genes within chromosomes, using 10,000 shifts and seed 20261013. Dog
tests cover the combined and significant sets; Human tests cover the complete
Human set. Empirical one-sided p-values use
`(1 + count(null mean ≤ observed mean)) / 10,001` and are unadjusted. Dog shifts
use batched indices with the original two-set RNG draw order. Cached shifts
are reused only when their ordered scores and membership, seed, number of
shifts, package versions and implementation checksum match.

The notebook regenerates trait-comparison rank CDFs, the combined six-trait
rank panels, enrichment tables, OCD threshold counts, and the additional
histograms and p-value/rank CDFs. Output filenames retain the existing names
inside this `toga/` directory. Earlier Ensembl-filtered 78-gene exports remain
outside it. The earlier executed notebook is preserved in
`../ensembl78_notebook.json`.

Section 5 is a separate constraint comparison and retains all 87 Dog and 48
Human ST6 genes. COL6A5 lacks a constraint estimate, leaving 86 Dog estimates.
Its methods and exports remain in [the constraint results](../../shet/METHODS.md).

## Validation and results

The notebook completed all 55 code cells without errors. All 20 fixture tests
passed, including TOGA selection and agreement of the batched circular shifts
with the original label-rolling algorithm. All 22 regenerated PNG figures were
visually inspected, along with the constraint figure.

`src/validate_toga_notebook.py` independently verifies all 36 Mann–Whitney tests,
genome-wide ranks and the six OCD Holm adjustments. It also checks that control
traits have no Holm values. All six traits' numeric window scores
match the established inputs; the selected 83 genes are covered in both windows
for every trait. The Human-set Mann–Whitney results and CDF source data match
the previous analysis. All 24 Dog circular-shift tests use 10,000 shifts and the
correct TOGA membership. See [the validation report](./toga_validation.json).

`src/validate_psychiatric_results.py` independently rescanned all four raw
psychiatric GWAS files, checking 664 Dog windows (83 genes × two windows × four
traits), genome-wide ranks, six Mann–Whitney tests per trait, all six OCD Holm
adjustments, and 80 sampled reverse liftover mappings. OCD threshold counts,
minima, lead SNPs and SNP counts were verified against the mapped SNP cache
inside the notebook. See [the raw-data validation report](./psychiatric_independent_validation.json).

The combined TOGA set remains enriched for OCD signal at both windows:

| Window | Genes | Mann–Whitney p | Holm p (six OCD tests) | Circular-shift p |
|---|---:|---:|---:|---:|
| ±1 Mb | 83 | 0.00106747 | 0.00320241 | 0.0390961 |
| ±500 kb | 83 | 4.01815×10⁻⁷ | 2.41089×10⁻⁶ | 0.000299970 |

Five of the six OCD tests have Holm p < 0.05. The significant and combined
sets pass at ±1 Mb, and all three sets pass at ±500 kb. The suggestive set at
±1 Mb has Holm p = 0.0825375. Control traits have unadjusted results only.
Full results are in
[mann_whitney_results.csv](./mann_whitney_results.csv) and
[circular_shift_results.csv](./circular_shift_results.csv).

This run rebuilt the BIP/ADHD/MDD/SCZ windows and confirmed their agreement with
the previous genome-wide scores. Subsequent default runs reuse verified caches.
The Python/package versions and input/code checksums are recorded in
`toga_validation.json`. This project has no Git repository, so reproducibility
is recorded through the saved source, executed notebook and file checksums.

The preceding TOGA outputs with 36-test Holm correction are preserved under
`holm36/`, including the executed notebook and methods. The optional
`ortholog_comparison_mann_whitney.csv` recomputes the historical 78-gene OCD
results with the same six-test family for comparison with the current 83-gene
results; it does not compare different correction families.

The OCD comparison remains exploratory because ST6 selection incorporated
the human OCD GWAS. Gene windows can overlap and share SNPs, and differences
in significance between traits do not establish a difference in enrichment.
