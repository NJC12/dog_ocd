"""Write plain-text explanations alongside the notebook's PDF figures."""

from pathlib import Path


TRAITS = {
    "CAD": "coronary artery disease",
    "BIP": "bipolar disorder",
    "ADHD": "attention-deficit/hyperactivity disorder",
    "MDD": "major depressive disorder",
    "SCZ": "schizophrenia",
}

WINDOWS = (
    "Each gene's association is its minimum GWAS SNP p-value within an inclusive "
    "window extending ±1 Mb or ±500 kb from its GENCODE transcription start site "
    "(TSS). Only covered autosomal protein-coding genes, with at least one valid "
    "SNP in the relevant window, receive an association score."
)
DOG_SET = (
    "The dog set contains 83 ST6 TOGA one-to-one human orthologs: 43 genes from "
    "significant dog GWAS loci and 40 from suggestive dog GWAS loci. The plotted "
    "associations are from human GWAS."
)
CAVEAT = (
    "These curves and distributions summarize gene-level association scores; "
    "they do not establish causality. Nearby genes can share SNPs and are not "
    "independent observations. Dog-gene selection used human OCD results, so "
    "the OCD comparison is exploratory rather than independent replication."
)


def write_plot_explanation(pdf_path):
    """Write a figure-specific explanation with the same stem as the PDF."""
    pdf_path = Path(pdf_path)
    stem = pdf_path.stem
    if stem == "shet_vs_ocd_500kb":
        paragraphs = [
            "Human constraint and OCD association in ST6 dog and human genes",
            "Each point represents one ST6 gene with an available human Agarwal "
            "posterior mean s_het and an OCD ±500 kb TSS-window association. "
            "This analysis uses the full ST6 groups, before the TOGA restriction: "
            "87 Dog and 48 Human genes. COL6A5 lacks a constraint estimate, "
            "leaving 86 Dog and 48 Human points. Both groups use human constraint estimates.",
            "The x-axis is −log10 of the minimum OCD SNP p-value within ±500 kb "
            "of the TSS; farther right indicates stronger OCD association. The "
            "y-axis is log10 of the posterior mean s_het; higher values indicate "
            "stronger estimated selection against heterozygous loss-of-function variants.",
            "Dark-blue and light-blue points denote significant and suggestive "
            "dog GWAS genes; orange points denote ST6 Human genes. Blue and orange "
            "squares mark the Dog and Human group medians, computed in the original "
            "units before transformation. The blue median combines both dog categories.",
            "The gray line is a pooled ordinary least-squares fit with an intercept "
            "across all 134 individual genes in the displayed coordinates. Its slope "
            "is printed on the plot; median squares are not additional observations.",
            "The top and right panels show Gaussian kernel density estimates for "
            "all Dog genes combined (blue) and Human genes (orange), using Scott's "
            "bandwidth rule separately for each group. Each gene has equal weight. "
            "These describe gene-level point estimates, not pooled posterior densities. "
            "Posterior uncertainty is not propagated into the scatter plot or fit.",
        ]
    else:
        panels = stem.endswith("_gene_rank_panels_all_traits")
        if panels:
            human = stem.startswith("human_")
            group = "human" if human else "dog"
            paragraphs = [
                f"Cumulative ranks of {group} genes across six human GWAS traits",
                "The left panel uses ±500 kb TSS windows and the right panel uses "
                "±1 Mb windows. Colors distinguish OCD, CAD, BIP, ADHD, MDD and "
                "SCZ; MDD uses a dashed line and SCZ a dash-dot line. CAD is coronary "
                "artery disease, BIP bipolar disorder, ADHD attention-deficit/"
                "hyperactivity disorder, MDD major depressive disorder, and SCZ schizophrenia.",
                "The Human set contains all 48 genes labeled Human in ST6, mapped "
                "by exact symbol to GENCODE." if human else DOG_SET,
            ]
        else:
            trait = next((key for key in TRAITS if stem.endswith("_" + key.lower())), "CAD")
            paragraphs = [
                f"Dog-gene associations in human OCD and {TRAITS[trait]} ({trait}) GWAS",
                DOG_SET,
            ]
        paragraphs.append(WINDOWS)
        if "rank" in stem:
            paragraphs.append(
                "The x-axis ranks each gene's minimum window p-value among all "
                "covered protein-coding genes for that trait and window. Rank 1 "
                "is the strongest association; ties receive average ranks. The "
                "y-axis is the cumulative fraction of the selected gene set at "
                "or below each rank. Curves that rise earlier indicate stronger "
                "associations within that trait's genome-wide rank distribution."
            )
            if "overlay" in stem:
                paragraphs.append(
                    "The four step curves show all 83 dog genes for OCD and the "
                    "comparison trait at both window sizes. The legend identifies "
                    "each trait/window combination. The black diagonal is a "
                    "uniform-rank reference."
                )
            elif not panels:
                paragraphs.append(
                    "Rows show OCD (top) and the comparison trait (bottom); columns "
                    "show ±1 Mb (left) and ±500 kb (right). Gray denotes all 83 dog "
                    "genes, orange the 43 significant genes, and blue the 40 suggestive genes."
                )
        elif stem.startswith("dog_gene_p_cdfs"):
            paragraphs.append(
                "Rows show OCD (top) and the comparison trait (bottom); columns "
                "show ±1 Mb (left) and ±500 kb (right). The x-axis is log10 of a "
                "p-value cutoff. The y-axis is the fraction of dog genes whose "
                "minimum window p-value is no greater than that cutoff. Curves "
                "that rise farther left indicate stronger associations. Gray "
                "denotes all 83 dog genes, orange the 43 significant genes, and "
                "blue the 40 suggestive genes. These are raw window minimum "
                "p-values, without adjustment for the number of SNPs per window."
            )
        elif stem.startswith("ocd_p_distributions"):
            paragraphs.append(
                "Rows show OCD ±1 Mb, OCD ±500 kb, comparison-trait ±1 Mb, and "
                "comparison-trait ±500 kb, from top to bottom. Columns show all "
                "covered protein-coding genes (gray), significant dog genes "
                "(orange), and suggestive dog genes (blue). Panel labels give gene counts. "
                "The gray group includes dog genes."
            )
            paragraphs.append(
                "The x-axis is −log10 of each gene's minimum window p-value; "
                "larger values indicate stronger association. Values above 15 "
                "are pooled in the final bin, labeled ≥15. The y-axis is gene "
                "density: each panel's histogram has unit area. The two window "
                "sizes share a density scale within each trait. Window minimum "
                "p-values are not adjusted for the number of SNPs per window."
            )
        else:
            raise ValueError(f"No plot explanation defined for {pdf_path.name}")
        paragraphs.append(CAVEAT)
    pdf_path.with_suffix(".txt").write_text("\n\n".join(paragraphs) + "\n", encoding="utf-8")
