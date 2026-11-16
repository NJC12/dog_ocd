# Removal candidates awaiting your decision

Nothing in this list has been deleted. Numbers identify groups you can approve
individually. Original inputs and historical analyses require a deliberate choice;
“not used by the current notebook” does not by itself make a file unnecessary.

1. **Local caches — safe to regenerate.** Root `.DS_Store` (Finder metadata) and all `*.pyc` files under
   `src/__pycache__/` and `tests/__pycache__/`. Python recreates these automatically;
   they are already excluded by `.gitignore`. They may reappear after tests run.

2. **Completed one-time Holm-policy migration — no longer needed for running the
   analysis.** It has already been applied. It edits fixed notebook cell positions
   and clears execution outputs, so it should not be rerun as part of the pipeline.

   [src/update_ocd_holm_notebook.py](../src/update_ocd_holm_notebook.py)

3. **Unused Agarwal gene/transcript map — optional external reference.** No current
   source script reads this file; the constraint analysis uses the shared Agarwal
   summary and posterior files instead. Removal would discard an original reference,
   not a generated cache. Include its `.txt.txt` provenance companion if removing it.

   [data/downloaded/constraint/gene_to_ensembl_id_map.txt](../data/downloaded/constraint/gene_to_ensembl_id_map.txt)

4. **OCD cohort-metadata workbook — optional source documentation.** It is not SNP
   summary statistics and no analysis script reads it. Retaining it preserves the
   publisher's study/cohort details. Include its `.xls.txt` companion if removing it.

   [basic.OCD_full_wo23andMe_190522.num.xls](../data/downloaded/metadata/basic.OCD_full_wo23andMe_190522.num.xls)

5. **Three legacy plotting scripts — optional if keeping only current TOGA plots.**
   Their descriptive analyses are included in the main notebook, but these scripts
   reproduce the historical 78-gene Ensembl-confirmed plots. Keep them if that
   reproduction workflow remains useful; the current notebook uses 83 TOGA genes.

   [src/plot_dog_gene_cdfs.py](../src/plot_dog_gene_cdfs.py)

   [src/plot_dog_gene_rank_cdfs.py](../src/plot_dog_gene_rank_cdfs.py)

   [src/plot_dog_gene_rank_overlay.py](../src/plot_dog_gene_rank_overlay.py)

6. **Standalone Human-set plotting script — optional independent check.** The
   notebook now contains the Human-set panels, while this script independently
   checks earlier saved ranks and audits that the notebook was not modified.
   Removing it would retain the notebook figure but lose that standalone check.

   [src/plot_human_gene_rank_panels.py](../src/plot_human_gene_rank_panels.py)

7. **Earlier pilot validation reports — optional historical records.** These record
   partial development-stage checks and are not read by the current pipeline.
   Later full-file validation reports exist. Keep them if the development history
   is useful; removing them would not remove the underlying analysis inputs.

   [psychiatric_pilot_validation.json](../notebook/results/ipynb/psychiatric_pilot_validation.json)

   [mdd_scz_pilot_validation.json](../notebook/results/ipynb/mdd_scz_pilot_validation.json)

## Files to retain

Do not remove `notebook/results/` wholesale. The default notebook and validation
tools read established per-gene scores and historical rank/Human-set tables from
it. The `mapped_grch38/` Parquet files are also used in default SNP-level threshold
checks. Rebuilding these would require the large GWAS files and configured chain.

The Ensembl REST/BioMart responses are downloaded legacy-audit evidence, not
disposable Python caches. The archived Ensembl78 and Holm36 notebooks preserve
different membership/correction policies. The burden analysis is absent from the
main notebook but has unique results. None is proposed for deletion here.

The PDF/TSV/TXT files in `plots_and_tables/` intentionally fulfill the requested
delivery exports. Full TSVs and displayed-table TSVs can overlap in content while
serving different purposes; both are retained.
