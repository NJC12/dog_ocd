# Repository Guidelines

## Project Structure & Module Organization

Analysis scripts are in `src/` (see `src/README.md` for notebook inclusion) and fixture tests in `tests/`. External inputs are in `data/downloaded/`: GWAS files in `gwas/`, GENCODE in `annotation/`, and the supplementary workbook in `supplementary/`. Each external file has a provenance/contents `.txt` companion. Transformed SNP caches and locally extracted metadata are in `data/generated/`. `notebook/results/` contains generated scientific results and historical scoring inputs; current notebook results go to `notebook/results/ipynb/toga/`. `dog_ocd_analysis.ipynb` stays at the root and also exports PDF/TSV files with matching explanations to `plots_and_tables/`. Keep generated files out of source and downloaded-input directories. See `data/README.md` for shared external input dependencies.

## Build, Test, and Development Commands

See README.md for notebook environment setup and `docs/analysis-guide.md` for script execution order. Install notebook dependencies from `config/requirements-notebook.txt` into a project-local environment. Run the existing tests from the project root with `.venv/bin/python -m unittest discover -s tests -v`. The notebook defaults to saved gene-window inputs and checks its recalculated ranks and Mann–Whitney results against the established CSVs.

## Coding Style & Naming Conventions

Use the standard formatter and linter for each language adopted by the project; commit their configuration alongside the code. Keep names descriptive and consistent within a language. For example, a Python module that checks inputs could be `src/validate_inputs.py`, with tests in `tests/test_validate_inputs.py`. Avoid unexplained abbreviations and unrelated formatting changes.

## Testing Guidelines

Tests use Python's `unittest` framework. Add focused tests with new behavior and bug fixes, using names that describe the behavior being checked. Run the full suite before opening a pull request and report the command and result in the request. For notebook changes, also execute all default cells and inspect generated figures.

## Commit & Pull Request Guidelines

This directory has no Git history from which to infer a commit convention. Use short, imperative commit subjects, such as `Add input validation`. Pull requests should explain the change, link a relevant issue when one exists, list verification steps, and include screenshots for visible interface changes.

## Terminology

- **Dog set:** in the executable notebook, the 83 ST6 TOGA one-to-one human orthologs: 43 dog-GWAS-significant and 40 dog-GWAS-suggestive genes. Legacy command-line results use 78 Ensembl-confirmed pairs.
- **Dog-set Holm family:** the six OCD tests (two windows × three dog sets). The five control traits report raw p-values without Holm adjustment; the Human-set family remains separate.
- **TSS window:** the inclusive interval extending ±1 Mb or ±500 kb from the GENCODE gene-level transcription start site.
- **Covered gene:** a gene with at least one valid GWAS SNP in the specified trait and TSS window.
- **Gene rank:** the genome-wide rank of a covered gene's minimum window p-value, with rank 1 strongest and average ranks for ties.
- **ST6 groups:** all 87 Dog and 48 Human genes in the supplementary table, before the dog-set orthology restriction; both use human constraint estimates.
- **Human set:** all 48 genes labeled Human in ST6, mapped by exact symbol to GENCODE; enrichment uses all other covered genes as background and a separate 12-test Holm family.
- **s_het rank:** pooled rank of available ST6 genes by human posterior mean s_het, with rank 1 highest and average ranks for ties.
