# Analysis source guide

The main analysis is the root-level [dog_ocd_analysis.ipynb](../dog_ocd_analysis.ipynb).
Run it from top to bottom using the project environment. Its current dog set is
the **83 ST6 TOGA one-to-one human orthologs**: 43 significant and 40 suggestive
dog GWAS genes. The Human set contains all 48 ST6 Human genes. The constraint
analysis uses the original 87 Dog and 48 Human ST6 groups, before TOGA filtering.

The notebook imports some modules and contains its own versions of other
functions. “Included” below describes whether the scientific analysis or utility
is represented in the notebook; it does not mean the standalone script is run.
Legacy scripts use the earlier **78 Ensembl-confirmed dog genes** (39 significant
and 39 suggestive) and can produce different gene-set results.

## Modules used directly by the notebook

| File | Responsibility | Inclusion in the notebook |
| --- | --- | --- |
| [compare_psychiatric.py](./compare_psychiatric.py) | Defines BIP, ADHD, MDD and SCZ input specifications; reads and validates GWAS chunks, lifts GRCh37 SNPs to GRCh38, scores TSS windows, and verifies or rebuilds score caches. | **Direct import.** Input previews and psychiatric window loading use this module. It also imports coordinate helpers from `compare_ocd.py`. |
| [toga_gene_set.py](./toga_gene_set.py) | Selects ST6 TOGA one-to-one dog genes, maps FAM155A to NALF1, assigns current membership without changing window scores, performs dog circular shifts, and applies Holm correction only to the six OCD dog-set tests. | **Direct import.** This defines the current 83-gene dog set and its multiple-testing policy. |
| [human_gene_enrichment.py](./human_gene_enrichment.py) | Maps the 48 ST6 Human genes, tests them against all other covered genes, applies a separate 12-test Holm family, and performs their circular-shift tests. | **Direct import.** Used in the Human-set analysis. |
| [compare_shet.py](./compare_shet.py) | Integrates Agarwal human s_het posteriors, ranks the complete ST6 groups, runs a two-sided Mann–Whitney comparison, and creates the constraint/OCD scatter plot, pooled regression and marginal densities. Writes results under `notebook/results/shet/`. | **Direct import and execution of `main()`.** Section 5 runs this analysis. The available estimates cover 86 Dog and 48 Human genes; COL6A5 is missing. |
| [notebook_exports.py](./notebook_exports.py) | Exports raw DataFrame/Series/Styler values as TSVs and writes table explanations, including column descriptions. | **Direct import.** Supports displayed tables and full analysis tables in `plots_and_tables/`. |
| [plot_explanations.py](./plot_explanations.py) | Writes figure-specific descriptions of gene sets, panels, axes, reference lines and interpretation alongside PDF figures. | **Direct import.** Each notebook PDF in `plots_and_tables/` receives a matching TXT file. |

The circular-shift modules use chromosome/SNP-count-stratified percentile ranks
for their test statistic. These differ from the genome-wide gene ranks shown in
the CDFs and Mann–Whitney tables.

## Standalone analyses also represented in the notebook

| File | Responsibility | Inclusion in the notebook |
| --- | --- | --- |
| [compare_ocd.py](./compare_ocd.py) | Provides chain parsing, coordinate lifting, GENCODE gene loading, OCD SNP validation, window minima, stratified ranks/circular shifts and distribution plots. Its standalone `main()` selects the legacy Ensembl-confirmed dog set and writes the established OCD results. | **Included through notebook implementations of the core helpers and scoring workflow.** Default notebook execution reads saved OCD windows; rebuilding is optional. Current dog membership comes from TOGA, and current circular shifts come from `toga_gene_set.py`. Coordinate helpers are also imported indirectly by `compare_psychiatric.py`. |
| [compare_cad.py](./compare_cad.py) | Verifies the CAD checksum, scans chromosome-sorted GRCh38 summary statistics, scores both window widths, and creates legacy CAD enrichment results and combined OCD/CAD distributions. | **Included through notebook implementations of checksum verification and CAD scoring.** Default execution reads saved CAD windows; the standalone `main()` is not called. |
| [preprocess_gwas.py](./preprocess_gwas.py) | Validates and lifts the full OCD GWAS, writing per-chromosome mapped SNP Parquet files and a liftover audit. | **Included through the notebook's `load_mapped_ocd()` workflow.** Cache rebuilding is optional; the standalone script is not called. |
| [rank_sum_enrichment.py](./rank_sum_enrichment.py) | Creates genome-wide covered-gene association ranks and one-sided Mann–Whitney tests for three dog sets. Its standalone `main()` processes legacy OCD/CAD inputs and adjusts all 12 tests together. | **Included through notebook `rank_and_test()` logic, with a different current correction policy.** The notebook evaluates six traits and applies Holm only to six OCD tests; control traits retain raw p-values. |
| [count_ocd_threshold_genes.py](./count_ocd_threshold_genes.py) | Checks saved OCD dog-gene window minima against mapped SNPs and counts genes below p < 5e-8 or p < 5e-6. The standalone script requires 78 legacy dog genes. | **Included through notebook implementations of classification and counting.** The notebook uses the 83-gene TOGA set and both TSS windows. |
| [plot_dog_gene_cdfs.py](./plot_dog_gene_cdfs.py) | Creates four-panel OCD/CAD minimum-p-value CDFs for the legacy combined, significant and suggestive dog sets. | **Included through notebook `plot_p_cdfs()`.** The notebook uses current TOGA membership and adds the psychiatric-trait comparisons. |
| [plot_dog_gene_rank_cdfs.py](./plot_dog_gene_rank_cdfs.py) | Creates four-panel OCD/CAD genome-wide gene-rank CDFs, requiring the 78/39/39 legacy gene counts. | **Included through notebook `plot_rank_cdfs()`.** The notebook uses 83/43/40 genes and adds the psychiatric-trait comparisons. |
| [plot_dog_gene_rank_overlay.py](./plot_dog_gene_rank_overlay.py) | Overlays four OCD/CAD trait-by-window CDF curves for all 78 legacy dog genes with a uniform-rank reference. | **Included through notebook `plot_rank_overlay()`.** The notebook uses 83 dog genes and produces separate OCD comparisons with CAD, BIP, ADHD, MDD and SCZ. |
| [plot_human_gene_rank_panels.py](./plot_human_gene_rank_panels.py) | Independently maps all 48 Human genes, checks saved six-trait genome-wide ranks, creates their two-window CDF panels and source data, and audits that the notebook was not modified. Writes to `notebook/results/human_gene_cdfs/`. | **The descriptive analysis is included through notebook `plot_all_trait_rank_panels()` with Human-set membership.** This standalone script is not called; it reads the earlier saved rank export and adds its own independent rank/audit checks. |

## Separate analysis and historical audit

| File | Responsibility | Inclusion in the notebook |
| --- | --- | --- |
| [compare_burden.py](./compare_burden.py) | Joins gene-level OCD pLoF SKAT-O burden results to legacy gene annotations, tests combined/significant/suggestive dog sets, adjusts the three one-sided tests with Holm, and creates rank/p-value CDFs and missing-gene audits. Writes to `notebook/results/burden/`. | **Not included.** This is a separate rare-variant gene-burden analysis, rather than the notebook's common-variant TSS-window analysis. |
| [audit_orthology.py](./audit_orthology.py) | Compares ST5/ST6 gene membership, maps ST6 symbols to GENCODE, queries/caches Ensembl BioMart and REST orthology, and records mapping decisions. The legacy pipeline uses its confirmed one-to-one decisions. | **The Ensembl audit is not included.** The notebook selects dog genes using the supplied ST6 TOGA labels and does not apply the Ensembl-confirmation filter. Existing legacy results remain relevant as score inputs and historical comparisons. |

## Execution, validation and migration tools

| File | Responsibility | Inclusion in the notebook |
| --- | --- | --- |
| [execute_analysis_notebook.py](./execute_analysis_notebook.py) | Runs all fixture tests, parses code cells, then executes and saves the root notebook. `--pilot` executes the first 28 cells and writes a separate pilot notebook. | **External runner.** Executes the notebook rather than contributing an additional analysis. |
| [validate_toga_notebook.py](./validate_toga_notebook.py) | Checks execution status, exact TOGA membership, scores against earlier inputs, genome-wide ranks, all 36 dog Mann–Whitney tests and current Holm policy, shift counts and unchanged Human-set outputs. Writes a validation JSON and temporary figure contact sheet. | **External validation.** Some assertions overlap notebook checks, but this validation program is not executed by the notebook. |
| [validate_psychiatric_results.py](./validate_psychiatric_results.py) | Rescans BIP/ADHD/MDD/SCZ files and independently checks dog-window minima/counts using direct interval masks, then verifies ranks, tests, score checksums and reverse-liftover samples. Supports `--output` for the result directory. | **External validation.** The independent raw-file rescan is not part of default notebook execution. |
| [compare_ortholog_results.py](./compare_ortholog_results.py) | Joins saved 78-gene and 83-gene results, harmonizes the old Mann–Whitney correction policy, and reports changes in effect sizes, p-values, shifts and restored genes. Writes comparison CSVs. | **Not included.** Historical sensitivity comparison, run separately after both result sets exist. |
| [update_ocd_holm_notebook.py](./update_ocd_holm_notebook.py) | Archives the previous 36-test correction outputs and edits notebook cells to apply the six-OCD-test policy. Uses fixed cell positions and clears execution outputs. | **One-time maintenance tool, not an analysis cell.** Its intended policy is already present in the notebook. Review before rerunning, especially after notebook edits. |

## Recommended execution

1. Install the notebook dependencies in the project-local environment using
   `config/requirements-notebook.txt`. Keep downloaded inputs and generated
   caches separate under `data/downloaded/` and `data/generated/`.
2. For the current analysis, run
   `.venv/bin/python src/execute_analysis_notebook.py`, or execute all cells in
   the root notebook. Default settings reuse established gene-window inputs and
   verified shift caches. Full notebook results go to
   `notebook/results/ipynb/toga/`; PDF/TSV delivery copies and explanations go to
   `plots_and_tables/`. Section 5 also writes its constraint results to
   `notebook/results/shet/`.
3. Run `.venv/bin/python src/validate_toga_notebook.py` after execution. Use
   `.venv/bin/python src/validate_psychiatric_results.py` when the independent
   full-file psychiatric validation is needed; it scans large source files.
4. Run `compare_burden.py` separately when reproducing the rare-variant analysis.
   Its input membership comes from the legacy saved OCD annotation, not the
   notebook's TOGA relabeling.

The historical OCD/CAD pipeline order is `audit_orthology.py` → optional
`preprocess_gwas.py` → `compare_ocd.py` → `compare_cad.py` →
`rank_sum_enrichment.py`, followed by the threshold and legacy plotting scripts.
It rewrites established legacy result files used by notebook comparisons; run
it deliberately when reproducing that analysis. `compare_psychiatric.py`,
`toga_gene_set.py`, `human_gene_enrichment.py`, `notebook_exports.py` and
`plot_explanations.py` are library modules rather than standalone analysis CLIs.
