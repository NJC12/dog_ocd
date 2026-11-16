# Dog OCD GWAS genes and human GWAS signal

This project compares human OCD, coronary artery disease (CAD), bipolar disorder (BIP), ADHD, major depressive disorder (MDD), and schizophrenia (SCZ) GWAS signal near human orthologs of dog compulsive-disorder GWAS genes with the autosomal protein-coding background. The OCD comparison is exploratory because the supplied ST6 belongs to a dog–human network analysis that used the same human OCD GWAS.

## Inputs

- `data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx`: ST6 dog genes and dog significance labels; ST5 lists their dog GWAS regions.
- `data/downloaded/gwas/daner_OCD_full_wo23andMe_190522.gz`: [Figshare file 53357735](https://figshare.com/articles/dataset/ocd2024/28707155?file=53357735); expected MD5 `9ab61cb7412c42a0ba207856af6da0f8`.
- `data/downloaded/annotation/gencode.v50.annotation.gtf.gz`: [GENCODE v50 comprehensive CHR GTF](https://www.gencodegenes.org/human/), GRCh38.
- `data/downloaded/gwas/GCST90132314.h.tsv.gz`: [GWAS Catalog harmonised CAD summary statistics](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz); expected MD5 `82e19c1817c52d5a94c8f42d9f48ccd0`. Its [metadata](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz-meta.yaml) specifies 1-based GRCh38 coordinates; a local metadata copy is in `data/downloaded/metadata/`.
- `data/downloaded/gwas/bip2024_multianc_no23andMe.gz`: [PGC bip2024, Figshare file 49760775](https://figshare.com/articles/dataset/bip2024/27216117?file=49760775), multi-ancestry meta-analysis excluding 23andMe; expected MD5 `6615ee8a7c6443df124848d17ba1336e`.
- `data/downloaded/gwas/ADHD2022_iPSYCH_deCODE_PGC.meta.gz`: [PGC adhd2022](https://figshare.com/articles/dataset/adhd2022/22564390), iPSYCH + deCODE + PGC meta-analysis (Demontis et al., *Nature Genetics*, 2023; [doi:10.1038/s41588-022-01285-8](https://doi.org/10.1038/s41588-022-01285-8)); expected MD5 `ea57120c720e880e34f67f57483dcbb4`.
- `data/downloaded/gwas/pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv.gz`: [PGC MDD2025, Figshare file 51487010](https://figshare.com/articles/dataset/GWAS_summary_statistics_for_major_depression_PGC_MDD2025_/27061255?file=51487010), multi-ancestry meta-analysis excluding 23andMe ([Adams et al., *Cell*, 2025](https://doi.org/10.1016/j.cell.2024.12.002)); expected MD5 `919829e44688fe7b0d745199a7893746`.
- `data/downloaded/gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz`: [PGC scz2022, Figshare file 34517861](https://figshare.com/articles/dataset/scz2022/19426775?file=34517861), primary multi-ancestry autosomal meta-analysis ([Trubetskoy et al., *Nature*, 2022](https://doi.org/10.1038/s41586-022-04434-5)); expected MD5 `e186621c9c290605ad7e19968bfc3dfc`.
- `~/data/liftover/hg19ToHg38.over.chain`: existing hg19-to-hg38 UCSC chain.

All six summary statistics files are in `data/downloaded/gwas/`; `data/generated/gwas_manifest.json` records paths, sizes and verified checksums. The BIP/ADHD study READMEs and MDD/SCZ embedded metadata specify hg19/GRCh37 coordinates. Downloaded metadata and Figshare snapshots are in `data/downloaded/metadata/`; extracted MDD/SCZ headers are in `data/generated/source_metadata/`. The notebook lifts these
SNPs to GRCh38 with the same chain as OCD, sorts the unordered input, and applies
the same inclusive TSS windows. `src/compare_psychiatric.py` implements their
validated reader and cache handling.

Downloaded inputs and orthology responses are retained in `data/downloaded/`; transformed SNP caches are in `data/generated/`. Derived tables, the figure, and the method summary are in `notebook/results/`.

## Run

### OCD pLoF burden enrichment

Run `.venv/bin/python src/compare_burden.py` to rank genes using the supplied
burden CSV's `P SKAT-O` column and repeat dog-gene Mann–Whitney enrichment tests.
Tables, rank/P-value CDFs, and validation details are saved in
`notebook/results/burden/`. The source CSV and these outputs are local-only and
are excluded from the public repository. The associations can be viewed at
https://allbyall.researchallofus.org/ using phenotype code
`obsessive_compulsive_disorder` (ancestry `meta`, pLoF, MAF < 0.001). See
[input publication notes](./publication-review.md) for download restrictions.

### Dog versus human gene constraint

Run `.venv/bin/python src/compare_shet.py` to compare Agarwal posterior mean
s_het for all 87 dog and 48 human genes in ST6. The script reads the posterior
density files in `~/data/constraint/shet/agarwal/` and the existing OCD gene
windows. Results, gene ranks, and the scatter plot are saved in
`notebook/results/shet/`; COL6A5 has no estimate and is excluded from the test.
Section 5 of the Jupyter notebook runs and displays this analysis.
See [methods and results](../notebook/results/shet/METHODS.md).

### Jupyter notebook

Open [`dog_ocd_analysis.ipynb`](../dog_ocd_analysis.ipynb) for the analysis organized
into libraries/basic functions, GWAS reading, window-based gene ranking, and
OCD-versus-CAD, OCD-versus-BIP, OCD-versus-ADHD, OCD-versus-MDD and OCD-versus-SCZ comparisons. Each comparison
includes Mann–Whitney tables, a single-panel four-curve rank CDF with its diagonal
reference, and the original chromosome/SNP-count-stratified circular-shift tests.
Dog-set selection follows the original ST6 TOGA labels: 83 one-to-one pairs,
comprising 43 dog-GWAS-significant and 40 dog-GWAS-suggestive genes. Ensembl
homology calls do not filter this notebook. `src/toga_gene_set.py` maps ST6
symbols to the human scoring annotation and supplies the circular shifts;
the BIP/ADHD/MDD/SCZ reader is imported from `src/compare_psychiatric.py`.
The notebook contains the other calculation and plotting functions.
See [TOGA methods and results](../notebook/results/ipynb/toga/METHODS.md).

Current notebook exports are saved in `notebook/results/ipynb/toga/`.
All plots generated by the notebook are also saved as PDFs in `plots and tables/`.
Displayed tables and full analysis tables are exported there as TSVs, retaining
numeric precision. Each PDF and TSV has a matching `.txt` explanation.
The earlier 78-gene Ensembl-filtered exports remain in `notebook/results/`
and directly under `notebook/results/ipynb/`; their original executed notebook
is preserved as `notebook/results/ipynb/ensembl78_notebook.json`.
Section 5 retains the separate constraint comparison of all 87 Dog and 48 Human
ST6 genes and exports to `notebook/results/shet/`.

Sections 4G and 4H add the MDD and SCZ comparisons. Section 4E shows all six trait CDFs in two panels (±500 kb and ±1 Mb) and collects
results in separate Mann–Whitney and circular-shift tables, one per window for
each test. The figure is exported as `dog_gene_rank_panels_all_traits.png` and
`.pdf`; the tables are `all_traits_mann_whitney_500kb.csv`,
`all_traits_mann_whitney_1mb.csv`, `all_traits_circular_shifts_500kb.csv` and
`all_traits_circular_shifts_1mb.csv`, with matching HTML versions, all under
`notebook/results/ipynb/toga/`.

Section 4F, **Combined CDF panels and enrichment tables, human genes**, repeats
the section 4E analysis for all 48 ST6 Human genes. Each window has one
Mann–Whitney and one circular-shift table with six rows. The Human-set background
is all other covered genes, and Holm adjustment spans a separate family of
12 tests (six traits × two windows). Circular shifts use the original
chromosome/SNP-count stratification, 10,000 shifts and seed 20261013; their
empirical p-values are unadjusted. `src/human_gene_enrichment.py` supplies these
tests. Shift caches include gene membership, ordered scores, seed, shift count
and implementation checksums, and respect `RUN_PERMUTATIONS` and
`REBUILD_WINDOWS`. Exports in `notebook/results/ipynb/toga/` include
`human_gene_rank_panels_all_traits.{png,pdf}`, `human_mann_whitney_results.csv`,
`human_circular_shift_results.csv`, and per-window CSV/HTML tables named
`all_traits_human_mann_whitney_{500kb,1mb}` and
`all_traits_human_circular_shifts_{500kb,1mb}`. These Human-set tests remain
exploratory because ST6 selection used the human OCD GWAS.

To prepare a project-local Python environment:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r config/requirements-notebook.txt
```

Select `.venv/bin/python` as the kernel in your notebook editor and run all cells
from the project root. The notebook selects the TOGA gene set directly from the
original workbook and reuses existing genome-wide window scores, replacing
the dog membership labels. It recomputes all ranks and Mann–Whitney tests,
checks that the OCD/CAD genome-wide ranks match the saved scores, and verifies
OCD threshold counts against the mapped SNP cache. BIP/ADHD/MDD/SCZ score caches
record source, chain and score-file checksums. Circular shifts are recalculated
for the TOGA sets; subsequent runs reuse results only when ordered scores,
gene membership, seed, shift count, package versions and implementation
checksums match. New exports go to `notebook/results/ipynb/toga/`.
Dog-set Holm adjustment in the notebook covers only the six OCD Mann–Whitney
tests (two windows × three dog sets). CAD, BIP, ADHD, MDD and SCZ are controls:
their raw p-values are reported, with blank Holm values in CSV and a dash in
displayed tables. The Human-set analysis retains its separate 12-test family.
The previous TOGA notebook/results with 36-test correction are archived under
`notebook/results/ipynb/toga/holm36/`. Original command-line OCD/CAD exports
retain their historical 12-test correction.

After running the notebook, `.venv/bin/python src/validate_psychiatric_results.py`
independently rescans the raw BIP/ADHD/MDD/SCZ files with direct interval masks to check
all 664 dog-gene windows (83 genes × two windows × four traits), genome-wide
ranks, Mann–Whitney tests, Holm adjustments
and sample reverse liftover mappings. The report is exported to
`notebook/results/ipynb/toga/psychiatric_independent_validation.json`.

Set `REBUILD_WINDOWS = True` to rescan all six GWAS inputs and rebuild gene scores;
this also reruns the original 10,000-shift permutation analysis. Set
`REBUILD_LIFTOVER = True` as well to regenerate the OCD Parquet cache using the
configured UCSC chain. `RUN_PERMUTATIONS = True` reruns permutations without
rebuilding windows. `SHOW_ADDITIONAL_PLOTS = True` displays the histograms and
four-panel CDFs for each OCD comparison and is enabled in the delivered notebook;
the rebuild and permutation-force options are off. The original command-line
scripts below retain the earlier Ensembl-filtered 78-gene analysis.

### Command-line scripts

From this directory with Python 3 and `numpy`, `pandas`, `openpyxl`, `pyarrow`, `matplotlib`, `requests`, `scipy`, and `statsmodels` installed:

```bash
mkdir -p data/downloaded/gwas data/downloaded/annotation
curl -L --fail -o data/downloaded/gwas/daner_OCD_full_wo23andMe_190522.gz https://ndownloader.figshare.com/files/53357735
curl -L --fail -o data/downloaded/annotation/gencode.v50.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz
md5 data/downloaded/gwas/daner_OCD_full_wo23andMe_190522.gz
python3 -m unittest discover -s tests -v
python3 src/preprocess_gwas.py
python3 src/audit_orthology.py
python3 src/compare_ocd.py
python3 src/compare_cad.py
python3 src/rank_sum_enrichment.py
python3 src/count_ocd_threshold_genes.py
python3 src/plot_dog_gene_cdfs.py
python3 src/plot_dog_gene_rank_cdfs.py
python3 src/plot_dog_gene_rank_overlay.py
```

`preprocess_gwas.py` validates the full OCD summary file and caches 1-based GRCh38 SNP coordinates in `data/generated/mapped_grch38/`. `audit_orthology.py` checks ST6 dog genes against Ensembl dog–human orthology and GENCODE gene IDs. `compare_ocd.py` calculates the OCD gene windows. `compare_cad.py` scans the GRCh38 CAD file by chromosome and applies the same ±1 Mb and ±500 kb TSS windows, five SNP-count strata per chromosome, and 10,000 chromosome-wise circular shifts. Its outputs are `notebook/results/per_gene_cad.csv` and `cad_summary.json`. The combined 12-panel figure is `notebook/results/ocd_p_distributions.png`. The random seed is fixed in `compare_ocd.py`.

`rank_sum_enrichment.py` ranks covered genes by their strongest window association and runs one-sided Mann–Whitney tests for dog significant, dog suggestive, and combined sets against the same non-dog covered genes. It writes `notebook/results/mann_whitney_gene_ranks.csv` and `mann_whitney_results.csv`. All 12 tests receive Holm-adjusted p-values; interpretation and limitations are in `notebook/results/METHODS.md`.

`count_ocd_threshold_genes.py` reports how many confirmed dog orthologs have at least one OCD SNP with P < 5×10⁻⁸ or P < 5×10⁻⁶ in either TSS window. It writes `notebook/results/ocd_threshold_counts.csv` and `ocd_threshold_dog_genes.csv` and verifies the per-gene values against the cached mapped SNPs.

`plot_dog_gene_cdfs.py` plots the cumulative fraction of dog genes whose strongest nearby association is at or below each p-value cutoff, for OCD and CAD at both window sizes. It writes `notebook/results/dog_gene_p_cdfs.png`.

`plot_dog_gene_rank_cdfs.py` plots the same dog-gene curves against each gene's rank among all covered protein-coding genes, with the strongest associations at the left. It writes `notebook/results/dog_gene_rank_cdfs.png`.

`plot_dog_gene_rank_overlay.py` puts the four trait-by-window rank curves for all 78 confirmed dog genes in one panel. It writes `notebook/results/dog_gene_rank_overlay.png`.
