# Human OCD and CAD GWAS signal around dog OCD GWAS genes

## Inputs and scope

The analysis uses the 6,791,632-row [no-23andMe OCD-3 summary file](https://figshare.com/articles/dataset/ocd2024/28707155?file=53357735), not the accompanying `.xls` cohort metadata. Its MD5 matched the publisher's `9ab61cb7412c42a0ba207856af6da0f8`. The [dataset readme](https://ndownloader.figshare.com/files/58046665) specifies GRCh37 SNP positions. Autosomal gene boundaries came from the [GENCODE v50 GRCh38 comprehensive annotation](https://www.gencodegenes.org/human/). The dog-gene labels came from ST6 of `NATURE_SupplementaryTables.xlsx`; ST5 and ST6 list the same 87 dog genes, of which ST6 marks 83 `1:1` after trimming whitespace and four `Not 1:1`.

The CAD analysis uses the supplied [GWAS Catalog GCST90132314 harmonised file](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz). Its MD5 matched the [published metadata](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz-meta.yaml): `82e19c1817c52d5a94c8f42d9f48ccd0`. The harmonised `chromosome` and `base_pair_location` columns are 1-based GRCh38, so no further liftover was applied. The original `markername` column encodes GRCh37 positions and was not used for interval searches. All 19,520,119 rows had valid autosomal chromosome, positive integer position, and p-value in [0,1]; none had p=0. Chromosome grouping and within-chromosome position order were checked during the scan. The metadata identifies this as the primary discovery meta-analysis of [Aragam et al. 2022](https://pubmed.ncbi.nlm.nih.gov/36474045/).

## Coordinate and orthology checks

SNP BP values were treated as 1-based. The existing `~/data/liftover/hg19ToHg38.over.chain` was applied to their corresponding 0-based positions, then converted back to 1-based GRCh38. Overlapping chain mappings were classified as ambiguous. A sample of 20 mapped SNPs reverse-mapped uniquely to its source coordinate. Of 6,791,632 input rows, all had valid autosomal chromosome, position, and p-value fields; 2,794 were unmapped, zero ambiguous, 771 mapped outside the autosomes, and 6,788,067 mapped to autosomes (including six that changed autosomal chromosome). The exact count audit is in `gwas_liftover_audit.json`.

The human ortholog of each dog gene was checked by stable Ensembl gene ID against Ensembl dog–human orthology records (bulk BioMart export with individual REST fallbacks). `FAM155A` was resolved through the `NALF1` stable human ID `ENSG00000204442`. Of the 83 ST6 `1:1` dog genes, 78 had a matching current Ensembl one-to-one ortholog, three had conflicting orthology (`EMC2`, `PPA2`, `SH3GL3`), and two could not be verified (`PLXNA2`, `PROSER3`). The four ST6 non-one-to-one genes were excluded. All decisions and Ensembl IDs are in `orthology_audit.csv` and `dog_ortholog_decisions.csv`. The primary dog set contains only the 78 confirmed pairs.

## Gene windows and comparison

For each autosomal GENCODE protein-coding gene, the single TSS is the annotated gene start on `+` strands and gene end on `−` strands. Both ±1 Mb and ±500 kb windows include both endpoints of `[max(1,TSS−width), TSS+width]`. At each distance and for each trait, the smallest valid SNP p-value, its rsID (or chromosome:position:alleles when rsID is missing) and GRCh38 position, and the number of covered SNP rows were recorded. Genes without SNP coverage remain missing. The per-gene tables are `per_gene_ocd.csv` and `per_gene_cad.csv`; ±500 kb columns end in `_500kb`.

The figure `ocd_p_distributions.png` shows unadjusted per-gene minimum p-values as `−log10(p)` histograms in 12 panels. Rows are OCD ±1 Mb, OCD ±500 kb, CAD ±1 Mb, and CAD ±500 kb. Columns show all covered genes, dog-GWAS-significant genes, and dog-GWAS-suggestive genes. Every panel shares the same x-axis from 0 to 15 and the same 0.5-unit bins; values at or above 15 are pooled in the final bin (1,218 CAD background genes at ±1 Mb and 694 at ±500 kb; no dog genes). Trait-specific vertical scales make the CAD panels legible; gene densities are comparable within a trait.

The figure `dog_gene_p_cdfs.png` shows empirical cumulative distributions of the same per-gene minimum p-values for the 78 confirmed dog orthologs. Its rows are OCD and CAD; its columns are ±1 Mb and ±500 kb. Each panel has curves for all 78 dog genes, the 39 dog-GWAS-significant genes, and the 39 dog-GWAS-suggestive genes. At x = log₁₀(P cutoff), the y-value is the fraction of genes in that set with a window minimum p-value ≤ 10ˣ. Both axes are shared across the four panels; the x-axis runs from −14 to 0. All dog genes have SNP coverage in both windows for both traits. Reproduce the plot with `python3 src/plot_dog_gene_cdfs.py`.

The alternate figure `dog_gene_rank_cdfs.png` uses the same four panels and dog sets, but its x-axis is each gene's rank among all covered autosomal protein-coding genes for that trait and window. Genes are ordered by their smallest window p-value, with rank 1 strongest and tied p-values receiving average ranks. At rank *r*, the y-value is the fraction of that dog set with rank ≤ *r*. The common x-axis extends to 19,166, the largest covered gene count. Ranks come from `mann_whitney_gene_ranks.csv`; reproduce the figure with `python3 src/plot_dog_gene_rank_cdfs.py`.

The single-panel figure `dog_gene_rank_overlay.png` shows only the combined set of 78 dog genes, with four colored curves for OCD and CAD at both window sizes. Each curve uses the rank of the gene among all covered protein-coding genes for its own trait and window. Thus the y-value at rank *r* is the fraction of the 78 dog genes with rank ≤ *r*. The four ranked backgrounds contain 19,119–19,166 genes. A straight black reference line connects the plot's minimum rank and fraction `(1, 0)` to its maximum rank and fraction `(19,166, 1)`, showing a uniform spread of ranks. Reproduce it with `python3 src/plot_dog_gene_rank_overlay.py`.

| Trait | Window | Covered background genes | Covered dog genes | Median dog minimum p | Median background minimum p | Dog genes with p < 5×10⁻⁸ |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| OCD | ±1 Mb | 19,144 | 78 | 1.64×10⁻⁴ | 3.50×10⁻⁴ | 0 |
| OCD | ±500 kb | 19,119 | 78 | 3.67×10⁻⁴ | 8.38×10⁻⁴ | 0 |
| CAD | ±1 Mb | 19,166 | 78 | 1.76×10⁻⁵ | 1.22×10⁻⁵ | 16 |
| CAD | ±500 kb | 19,166 | 78 | 1.09×10⁻⁴ | 4.76×10⁻⁵ | 9 |

### OCD threshold counts for dog genes

Among the 78 confirmed dog–human one-to-one orthologs, a gene counts if at least one mapped OCD SNP falls within the inclusive TSS window and has `P` strictly below the threshold. Equivalently, its minimum window p-value is below the threshold. The ±1 Mb and ±500 kb results were recalculated from the cached, lifted SNP rows and checked against every dog-gene minimum, lead SNP, and SNP count in `per_gene_ocd.csv`. These are gene counts, not independent GWAS loci; overlapping gene windows may count the same SNP more than once.

| TSS window | Dog set | P < 5×10⁻⁸ | Additional at P < 5×10⁻⁶ | Total at P < 5×10⁻⁶ |
| --- | --- | ---: | ---: | ---: |
| ±1 Mb | Significant | 0 | 3 | 3 |
| ±1 Mb | Suggestive | 0 | 2 | 2 |
| ±1 Mb | Combined | 0 | 5 | 5 |
| ±500 kb | Significant | 0 | 1 | 1 |
| ±500 kb | Suggestive | 0 | 1 | 1 |
| ±500 kb | Combined | 0 | 2 | 2 |

At ±1 Mb, the five genes below 5×10⁻⁶ are `RNF144A` (1.427×10⁻⁶), `RORB` (2.432×10⁻⁶), `NTRK2` (3.987×10⁻⁶), `CDH8` (4.280×10⁻⁶), and `SLIT2` (4.319×10⁻⁶). At ±500 kb, only `RORB` and `CDH8` remain. The complete per-gene audit and summary counts are `ocd_threshold_dog_genes.csv` and `ocd_threshold_counts.csv`; reproduce them with `python3 src/count_ocd_threshold_genes.py`.

For inference at each distance, each covered gene was assigned to one of five SNP-count strata within its chromosome, then ranked by minimum p-value within that chromosome and stratum. Smaller ranks indicate stronger signal. Ten thousand chromosome-wise circular shifts of dog labels over TSS-ordered covered genes supplied the null distribution. Each empirical one-sided p-value uses `(1 + count(null ≤ observed))/(10,001)`. The permutation seed was 20261013.

| Trait | Window | Dog set | Observed mean rank | Null mean rank | Empirical p |
| --- | --- | --- | ---: | ---: | ---: |
| OCD | ±1 Mb | All 78 confirmed pairs | 0.4126 | 0.5028 | 0.0247 |
| OCD | ±1 Mb | 39 dog-GWAS-significant pairs | 0.3757 | 0.5034 | 0.0478 |
| OCD | ±500 kb | All 78 confirmed pairs | 0.3512 | 0.5024 | 0.000300 |
| OCD | ±500 kb | 39 dog-GWAS-significant pairs | 0.3373 | 0.5038 | 0.0101 |
| CAD | ±1 Mb | All 78 confirmed pairs | 0.5325 | 0.5027 | 0.7192 |
| CAD | ±1 Mb | 39 dog-GWAS-significant pairs | 0.5630 | 0.5015 | 0.7584 |
| CAD | ±500 kb | All 78 confirmed pairs | 0.5905 | 0.5024 | 0.9683 |
| CAD | ±500 kb | 39 dog-GWAS-significant pairs | 0.6553 | 0.5025 | 0.9814 |

Full precision and coverage counts are in `summary.json` and `cad_summary.json`.

## Mann–Whitney rank-sum comparison

For the rank-sum analysis, “maximum p-value” was clarified to mean the strongest association: the largest `−log10(p)` (equivalently, the smallest numerical p-value) among SNPs in each TSS window. For each trait and distance, all covered autosomal protein-coding genes were ranked by this value, with rank 1 strongest and tied values assigned their average rank. Three dog sets were tested: the 39 dog-GWAS-significant genes, the 39 dog-GWAS-suggestive genes, and their union of 78 confirmed dog orthologs. **Every test uses the same other covered protein-coding genes as its reference group, excluding all 78 dog genes.** Uncovered genes are excluded. The one-sided asymptotic Mann–Whitney U test asks whether dog genes have larger `−log10(p)` values; the U statistic incorporates ties. `U/(n_dog × n_other)` estimates the probability that a randomly selected gene from the specified dog set has stronger signal than a randomly selected non-dog gene, counting ties as one half. Holm-adjusted p-values cover all 12 trait-by-window-by-set tests. Full gene ranks and results are in `mann_whitney_gene_ranks.csv` and `mann_whitney_results.csv`.

| Trait | Window | Dog set | Dog genes | Non-dog genes | U/(n dog × n non-dog) | One-sided p | Holm p |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| OCD | ±1 Mb | Significant | 39 | 19,066 | 0.6477 | 0.000705 | 0.00564 |
| OCD | ±1 Mb | Suggestive | 39 | 19,066 | 0.5683 | 0.0700 | 0.490 |
| OCD | ±1 Mb | Combined | 78 | 19,066 | 0.6080 | 0.000488 | 0.00439 |
| OCD | ±500 kb | Significant | 39 | 19,041 | 0.6682 | 0.000139 | 0.00153 |
| OCD | ±500 kb | Suggestive | 39 | 19,041 | 0.6540 | 0.000438 | 0.00438 |
| OCD | ±500 kb | Combined | 78 | 19,041 | 0.6611 | 4.37×10⁻⁷ | 5.25×10⁻⁶ |
| CAD | ±1 Mb | Significant | 39 | 19,088 | 0.4192 | 0.960 | 1.000 |
| CAD | ±1 Mb | Suggestive | 39 | 19,088 | 0.4850 | 0.627 | 1.000 |
| CAD | ±1 Mb | Combined | 78 | 19,088 | 0.4521 | 0.928 | 1.000 |
| CAD | ±500 kb | Significant | 39 | 19,088 | 0.3692 | 0.998 | 1.000 |
| CAD | ±500 kb | Suggestive | 39 | 19,088 | 0.4888 | 0.596 | 1.000 |
| CAD | ±500 kb | Combined | 78 | 19,088 | 0.4290 | 0.985 | 1.000 |

The dog-significant OCD set shows enrichment at both distances, while the dog-suggestive OCD set shows enrichment at ±500 kb only after correction. The combined-set results agree in direction with the earlier chromosome-wise circular-shift analysis: OCD shows dog-gene enrichment and CAD does not. The rank-sum test does not account for gene-window overlap, local LD, SNP-count differences, or dog-gene genomic clustering, so its p-values should not replace the permutation results as the primary inference. The OCD comparison remains exploratory because the dog-gene selection used human OCD results.

## Interpretation

The confirmed dog genes show enrichment for stronger nearby human OCD GWAS signal under both window definitions, with a smaller empirical p-value at ±500 kb. The same dog set does not show enrichment for nearby CAD signal under the corresponding one-sided tests; its observed mean ranks are above the CAD null means. CAD is a different and much larger GWAS, so the raw minimum p-values are not directly comparable with OCD across traits. Windows overlap and can share lead SNPs, so the per-gene values are correlated; the chromosome-wise shifts preserve dog-gene spacing in gene order but do not establish causal gene assignments or independence from local LD. The ±500 kb reanalysis and subsequent CAD comparison were requested after the original OCD result. The distances, traits, and significant-only sensitivity results have no multiple-testing correction.

This is **exploratory, not independent replication**. ST6 was assembled for a dog–human network analysis that incorporates the [human OCD GWAS](https://www.nature.com/articles/s41588-025-02189-z), and the dog-gene selection cannot be demonstrated to be independent of those human results from the supplied table. The [associated dog–human study description](https://data.europa.eu/data/datasets/https-doi-org-10-17044-scilifelab-33339309?locale=en) explicitly describes combining the dog and human GWAS for network colocalization.

## Reproduction

Run the commands in the project [README](../../README.md). The code is in `src/`, fixture tests are in `tests/`, and `data/mapped_grch38/` caches the lifted OCD summary file. Tests cover TSS strand, inclusive window boundaries, nested ±500 kb/±1 Mb windows, overlapping and empty windows, one-based liftover, CAD parsing with p=0, invalid rows, and missing gene coverage, and rank-sum grouping, ties, and exclusion of uncovered genes. All seven fixture tests passed. Full-table checks confirmed that every ±500 kb SNP count is no larger than its ±1 Mb count, and every covered ±500 kb minimum p-value is at least its ±1 Mb counterpart, for both traits.
