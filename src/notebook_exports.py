"""Plain-text companion files and TSV exports for the executable notebook."""

from pathlib import Path

import pandas as pd
from pandas.io.formats.style import Styler


COLUMN_DESCRIPTIONS = {
    "trait": "Human GWAS trait (OCD, CAD, BIP, ADHD, MDD or SCZ).",
    "window": "TSS window: 1mb is ±1 Mb; 500kb is ±500 kb, inclusive.",
    "dog_set": "Significant (43), suggestive (40), or combined (83) TOGA dog genes.",
    "gene_id": "GENCODE human Ensembl gene identifier.",
    "symbol": "Human gene symbol.",
    "chrom": "Autosomal chromosome.",
    "start": "GENCODE gene start, 1-based GRCh38.",
    "end": "GENCODE gene end, 1-based GRCh38.",
    "strand": "GENCODE gene strand (+ or −).",
    "tss": "Gene transcription start site, 1-based GRCh38; start on + strand, end on − strand.",
    "tss_grch38": "Gene transcription start site, 1-based GRCh38.",
    "min_p": "Minimum SNP p-value in the relevant window (±1 Mb if no window column).",
    "min_p_500kb": "Minimum SNP p-value in the ±500 kb TSS window.",
    "n_snps": "Number of valid SNPs in the relevant window (±1 Mb if no window column).",
    "n_snps_500kb": "Number of valid SNPs in the ±500 kb TSS window.",
    "lead_snp": "SNP attaining the minimum window p-value (±1 Mb if no window column).",
    "lead_snp_500kb": "SNP attaining the minimum ±500 kb window p-value.",
    "lead_pos_grch38": "Lead SNP position, 1-based GRCh38.",
    "lead_pos_grch38_500kb": "Lead SNP position for the ±500 kb window, 1-based GRCh38.",
    "rank": "Genome-wide rank of minimum window p-value; 1 is strongest, ties have average ranks.",
    "dog": "True if the gene belongs to the 83-gene TOGA dog set.",
    "dog_significant": "True if the gene belongs to the 43-gene significant dog subset.",
    "dog_suggestive": "True if the gene belongs to the 40-gene suggestive dog subset.",
    "n_dog": "Number of covered genes in the selected dog subset.",
    "n_human": "Number of covered ST6 Human genes.",
    "n_other": "Number of background genes: non-dog for dog tests, all other genes for Human-set tests.",
    "u_statistic": "Mann–Whitney U for the selected set versus its background, using association strength.",
    "probability_dog_stronger": "Probability a dog gene has stronger association than a non-dog gene, counting ties as one-half.",
    "probability_human_stronger": "Probability an ST6 Human gene has stronger association than a background gene, counting ties as one-half.",
    "one_sided_p": "One-sided asymptotic Mann–Whitney p-value testing stronger association in the selected set.",
    "holm_adjusted_p": "Holm-adjusted p-value; dog tests use six OCD tests only (control traits blank); Human-set tests use a separate 12-test family.",
    "holm_family_tests": "Number of tests in the Holm multiple-testing family.",
    "width_each_side_bp": "TSS window extension on each side, in base pairs.",
    "cdf_fraction": "Cumulative fraction of selected genes up to this row's association rank.",
    "gene_set": "Selected gene-set label.",
    "mean rank": "Observed mean chromosome/SNP-count-stratified signal-rank percentile; smaller means stronger association.",
    "null mean rank": "Mean selected-set signal-rank percentile across the circular-shift null replicates.",
    "empirical p": "One-sided circular-shift p-value with a plus-one correction, testing smaller mean rank.",
    "permutations": "Number of within-chromosome circular-shift replicates.",
    "shifts": "Number of within-chromosome circular-shift replicates.",
    "circular_shift_p": "One-sided empirical p-value from 10,000 within-chromosome circular shifts.",
    "genome_wide": "Count with minimum window p-value < 5e-8.",
    "additional_at_5e-6": "Count with 5e-8 ≤ minimum window p-value < 5e-6.",
    "neither": "Count with minimum window p-value ≥ 5e-6.",
    "total_below_5e-6": "Count with minimum window p-value < 5e-6, including genome-wide signals.",
    "genes_tested": "Number of dog genes evaluated.",
    "threshold_category": "Mutually exclusive genome_wide, additional_at_5e-6, or neither category.",
    "st6_symbol": "Gene symbol as recorded in supplementary table ST6.",
    "group": "ST6 Dog or Human gene-set label; both use human constraint estimates.",
    "shet_rank": "Pooled rank by posterior mean human s_het, highest first, average ranks for ties.",
    "shet_posterior_mean": "Mean of the human Agarwal s_het posterior, integrated on the linear s_het scale.",
    "log10_shet": "Base-10 logarithm of the posterior mean s_het.",
    "neg_log10_p_500kb": "Negative base-10 logarithm of the minimum OCD ±500 kb window p-value.",
    "shet_status": "Whether an Agarwal constraint estimate was available.",
    "dog_significance": "Significant or suggestive ST6 dog GWAS category.",
    "slope": "Pooled OLS slope in the displayed log-transformed coordinates.",
    "intercept": "Pooled OLS intercept in the displayed log-transformed coordinates.",
    "r_squared": "Fraction of displayed y-axis variance explained by the pooled OLS fit.",
    "human_symbol": "Mapped human ortholog gene symbol.",
    "human_gene_id": "Mapped GENCODE human Ensembl gene identifier.",
    "gene*": "Original ST6 dog gene symbol.",
    "dogcd significance": "ST6 significant or suggestive dog GWAS category.",
    "st6 toga label": "Orthology label supplied in ST6.",
    "genes": "Number of ST6 genes with this orthology label.",
    "max_neg_log10_p": "Negative log10 of the minimum SNP p-value in this trait/window.",
    "median_dog_min_p": "Median minimum-window p-value among the selected dog genes.",
    "median_human_min_p": "Median minimum-window p-value among the 48 ST6 Human genes.",
    "median_other_min_p": "Median minimum-window p-value among background genes.",
    "median_dog_rank": "Median genome-wide association rank among selected dog genes.",
    "median_human_rank": "Median genome-wide association rank among selected Human genes.",
    "median_other_rank": "Median genome-wide association rank among background genes.",
    "covered genes": "Number of protein-coding genes with at least one valid SNP in this window.",
    "covered dog genes": "Number of covered genes in the 83-gene TOGA dog set.",
    "index": "Row label retained from the source table, such as an audit category.",
    "snp rows": "Number of SNP records in the indicated audit category.",
    "toga": "Orthology label as supplied in ST6.",
    "agarwal_symbol": "Gene symbol used to retrieve the Agarwal posterior.",
    "shet_mapping_evidence": "Evidence for the ST6-to-Agarwal gene-symbol mapping.",
    "ensembl_transcript_id": "Ensembl transcript identifier supplied by Agarwal.",
    "agarwal_chromosome_group": "Chromosome category supplied by Agarwal.",
    "posterior_column": "Column identifying this gene in its source posterior-density file.",
    "density_integral": "Numerical integral of the supplied posterior density on its log10(s_het) grid.",
    "posterior_file": "Source file containing the gene's posterior density.",
    "in_previous_78_dog_set": "Legacy 78-gene Ensembl-confirmed dog membership; not the current TOGA selection.",
    "u_dog": "Mann–Whitney U for Dog versus Human posterior mean s_het.",
    "p_two_sided": "Two-sided asymptotic Mann–Whitney p-value comparing Dog and Human s_het.",
    "probability_dog_greater": "Probability a Dog gene's posterior mean s_het exceeds a Human gene's, counting ties as one-half.",
    "median_shet_dog": "Median posterior mean s_het among available ST6 Dog genes.",
    "median_shet_human": "Median posterior mean s_het among ST6 Human genes.",
    "median_rank_dog": "Median pooled s_het rank among available ST6 Dog genes.",
    "median_rank_human": "Median pooled s_het rank among ST6 Human genes.",
    "n": "Number of individual genes used in the fit.",
    "method": "Statistical method and settings.",
    "n_plotted": "Number of genes plotted with available constraint and OCD association estimates.",
    "missing_shet_genes": "ST6 genes lacking an Agarwal posterior.",
    "variable": "Displayed coordinate for the density estimate: neg_log10_p_500kb or log10_shet.",
    "value": "Evaluation point in the transformed coordinate named by variable.",
    "density": "Gaussian KDE density at value, with Scott's bandwidth and equal weight per gene within each group.",
}


def export_table(table, output_dir, stem, explanation, index=False):
    """Export unrounded table values and a contextual plain-text explanation."""
    if isinstance(table, Styler):
        table = table.data
    if isinstance(table, pd.Series):
        table = table.to_frame()
    table = table.copy()
    if index or table.index.name is not None or not pd.api.types.is_integer_dtype(table.index.dtype):
        table = table.reset_index()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{stem}.tsv"
    table.to_csv(path, sep="\t", index=False)
    descriptions = []
    aliases = {"dog genes": "n_dog", "human genes": "n_human", "non-dog genes": "n_other",
               "other genes": "n_other", "u": "u_statistic", "one-sided p": "one_sided_p",
               "holm p": "holm_adjusted_p", "p(dog stronger)": "probability_dog_stronger",
               "p(human stronger)": "probability_human_stronger"}
    for column in table.columns:
        key = str(column).lower()
        key = aliases.get(key, key)
        meaning = COLUMN_DESCRIPTIONS.get(key)
        if stem.startswith('shet') and key in ('n_dog', 'n_human'):
            meaning = "Number of ST6 genes in this group with available human constraint estimates."
        if meaning:
            descriptions.append(f"{column}: {meaning}")
        else:
            descriptions.append(f"{column}: {str(column).replace('_', ' ')}; see the table context above.")
    text = (
        explanation + "\n\n"
        f"Rows: {len(table):,}. Tab-separated UTF-8 text with one header row. "
        "Numeric values retain computational precision rather than notebook display rounding. "
        "Blank fields indicate missing or inapplicable values.\n\n"
        "Columns\n" + "\n".join(descriptions) + "\n"
    )
    path.with_suffix(".txt").write_text(text, encoding="utf-8")


def saved_table_explanation(stem):
    """Explain full analysis exports by their established notebook filenames."""
    if stem.startswith("per_gene_"):
        return (f"Per-gene human {stem.removeprefix('per_gene_').upper()} GWAS scores for all "
                "GENCODE autosomal protein-coding genes, including uncovered genes. Windows "
                "extend ±1 Mb and ±500 kb from the gene TSS, inclusive, using GRCh38 coordinates. "
                "Minimum SNP p-values are not adjusted for the number of SNPs per window. "
                "Dog labels identify the 83 TOGA one-to-one human orthologs (43 significant, 40 suggestive).")
    if "shet" in stem:
        return ("ST6 human constraint analysis: all 87 Dog and 48 Human genes, before TOGA "
                "restriction, using human Agarwal s_het estimates. COL6A5 lacks an estimate, "
                "leaving 86 Dog and 48 Human genes for analysis. Posterior means are integrated "
                "on the linear scale from densities supplied on a log10 grid. The group test "
                "is two-sided asymptotic Mann–Whitney with tie and continuity corrections. "
                "The regression pools individual genes in the plotted log-transformed coordinates. "
                "Marginal densities use Gaussian KDE with Scott's bandwidth separately by group; "
                "they describe point estimates and do not propagate posterior uncertainty.")
    if "circular" in stem:
        return ("Circular-shift enrichment tests for " +
                ("the 48 ST6 Human genes" if "human" in stem else "the significant and combined TOGA dog sets") +
                ". Each trait/window uses 10,000 independent within-chromosome circular shifts "
                "of gene-set membership among covered genes ordered by genomic position. The "
                "statistic is the mean signal-rank percentile within chromosome × SNP-count strata "
                "(smaller is stronger), rather than the genome-wide ranks used in CDFs. Empirical "
                "one-sided p-values use a plus-one correction. These p-values are unadjusted.")
    if "mann_whitney" in stem and "ranks" not in stem:
        return ("One-sided Mann–Whitney enrichment tests of stronger minimum-window association "
                "in the selected genes, using genome-wide ranks with average ties. " +
                ("The 48 ST6 Human genes are compared with all other covered genes. Holm correction "
                 "uses a separate 12-test family (six traits × two windows)." if "human" in stem else
                 "The 43 significant, 40 suggestive and 83 combined TOGA dog genes are each compared "
                 "with non-dog covered genes. Holm correction uses only the six OCD tests "
                 "(two windows × three dog sets); control traits report raw p-values.") +
                " Tests use asymptotic p-values with tie and continuity corrections. Nearby genes "
                "can share SNPs; OCD comparisons involving the dog set are exploratory because "
                "dog-gene selection used human OCD results.")
    if "threshold" in stem:
        return ("TOGA dog-gene OCD window association thresholds, independently checked against "
                "the mapped SNP cache. Each of 83 genes is evaluated at ±1 Mb and ±500 kb from "
                "its TSS. Categories are genome-wide (p < 5e-8), additional at 5e-6 "
                "(5e-8 ≤ p < 5e-6), and neither (p ≥ 5e-6). Combined counts sum significant "
                "and suggestive dog subsets. Minimum-window p-values are not adjusted for window SNP count.")
    if "gene_set" in stem:
        return ("Gene-set membership and annotation from supplementary table ST6. " +
                ("The Human set contains all 48 ST6 Human genes, mapped by exact symbol to GENCODE."
                 if "human" in stem else
                 "The dog set contains 83 TOGA one-to-one human orthologs: 43 significant "
                 "and 40 suggestive dog GWAS genes, selected from the original 87 ST6 Dog genes."))
    return ("Gene-level association ranks or cumulative-rank source data across human GWAS "
            "traits and inclusive ±1 Mb/±500 kb TSS windows. Covered genes are ranked genome-wide "
            "by minimum window p-value, strongest first, with average ranks for ties. The dog "
            "set uses 83 TOGA one-to-one human orthologs; the Human set uses all 48 ST6 Human genes.")
