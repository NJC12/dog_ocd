"""Plot the 48 ST6 Human genes using section 4E's six-trait CDF layout.

Run from the project root with the notebook environment. This reads saved gene
ranks and never writes to dog_ocd_analysis.ipynb or its established exports.
"""

import hashlib
import json
import platform
from pathlib import Path

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'notebook/results/human_gene_cdfs'
TRAITS = (
    ('OCD', 'OCD', '#d55e00'),
    ('CAD', 'CAD', '#0072b2'),
    ('BIP', 'Bipolar disorder (BIP)', '#cc79a7'),
    ('ADHD', 'ADHD', '#009e73'),
    ('MDD', 'Major depressive disorder (MDD)', '#6a3d9a'),
    ('SCZ', 'Schizophrenia (SCZ)', '#555555'),
)
PANELS = (('500kb', 'a  ±500 kb window'), ('1mb', 'b  ±1 Mb window'))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    notebook_path = ROOT / 'dog_ocd_analysis.ipynb'
    notebook_hash = sha256(notebook_path)
    workbook_path = ROOT / 'data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx'
    ranks_path = ROOT / 'notebook/results/ipynb/mann_whitney_gene_ranks.csv'

    st6 = pd.read_excel(workbook_path, sheet_name='ST6', header=1)
    human = st6.loc[st6.Species.eq('Human'), ['Gene*']].rename(
        columns={'Gene*': 'symbol'})
    if len(human) != 48 or human.symbol.isna().any() or human.symbol.duplicated().any():
        raise ValueError('Expected 48 unique Human genes in supplementary table ST6')
    symbols = set(human.symbol)
    ranked = pd.read_csv(ranks_path, dtype={'chrom': str})
    expected_groups = {(trait, window) for trait, _, _ in TRAITS for window, _ in PANELS}
    if set(zip(ranked.trait, ranked.window)) != expected_groups:
        raise ValueError('Saved ranks do not contain exactly six traits and two windows')

    selected_frames, audit = [], []
    human_ids = None
    frames = {}
    for trait, _, _ in TRAITS:
        for window, _ in PANELS:
            frame = ranked.loc[ranked.trait.eq(trait) & ranked.window.eq(window)].copy()
            if frame.empty or frame.gene_id.duplicated().any():
                raise ValueError(f'Missing or duplicate {trait} {window} gene ranks')
            if not frame['rank'].between(1, len(frame)).all():
                raise ValueError(f'Invalid {trait} {window} gene ranks')
            p = frame.min_p.to_numpy(dtype=float)
            if not np.isfinite(p).all() or not ((p > 0) & (p <= 1)).all():
                raise ValueError(f'Invalid {trait} {window} window p-values')
            if not frame.n_snps.gt(0).all():
                raise ValueError(f'Uncovered genes in saved {trait} {window} ranks')
            # Independently verify the saved ranks, including average tied ranks,
            # using sorted-array insertion positions rather than pandas.rank.
            scores = np.log10(p)
            sorted_scores = np.sort(scores)
            independently_ranked = (
                np.searchsorted(sorted_scores, scores, side='left')
                + np.searchsorted(sorted_scores, scores, side='right') + 1
            ) / 2
            np.testing.assert_array_equal(frame['rank'].to_numpy(), independently_ranked)

            selected = frame.loc[frame.symbol.isin(symbols)].copy()
            if len(selected) != 48 or selected.symbol.duplicated().any():
                raise ValueError(f'Expected 48 unique ST6 Human genes in {trait} {window}')
            if set(selected.symbol) != symbols:
                raise ValueError(f'ST6 Human symbols differ for {trait} {window}')
            ids = set(selected.gene_id)
            if human_ids is not None and ids != human_ids:
                raise ValueError('ST6 Human gene IDs differ across traits or windows')
            human_ids = ids
            selected = selected.sort_values(['rank', 'gene_id'], kind='mergesort')
            selected['cdf_fraction'] = np.arange(1, len(selected) + 1) / len(selected)
            selected['gene_set'] = 'ST6 Human'
            frames[trait, window] = selected
            selected_frames.append(selected)
            audit.append({'trait': trait, 'window': window,
                          'n_covered_background': len(frame), 'n_human': len(selected),
                          'median_human_rank': float(selected['rank'].median()),
                          'missing_human_symbols': []})

    max_rank = max(row['n_covered_background'] for row in audit)
    OUT.mkdir(parents=True, exist_ok=True)
    source_data = pd.concat(selected_frames, ignore_index=True)
    source_data = source_data[['trait', 'window', 'width_each_side_bp', 'gene_id',
                               'symbol', 'chrom', 'min_p', 'n_snps', 'rank',
                               'cdf_fraction', 'gene_set']]
    source_data.to_csv(OUT / 'human_gene_rank_cdf_source_data.csv', index=False)

    # Match section 4E; only the selected gene set and y-axis wording change.
    with plt.rc_context({'pdf.fonttype': 42}):
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.7), sharex=True, sharey=True)
        for ax, (window, title) in zip(axes, PANELS):
            for trait, label, color in TRAITS:
                selected = frames[trait, window]
                values = selected['rank'].to_numpy(dtype=float)
                fraction = selected.cdf_fraction.to_numpy(dtype=float)
                ax.step(np.r_[1, values, max_rank], np.r_[0, fraction, 1],
                        where='post', linewidth=2.2, color=color, label=label,
                        linestyle={'MDD': '--', 'SCZ': '-.'}.get(trait, '-'))
            ax.set_title(title, loc='left', fontsize=12, fontweight='bold')
            ax.set_xlim(1, max_rank)
            ax.set_ylim(0, 1.02)
            ax.set_xticks([1, 5_000, 10_000, 15_000, max_rank])
            ax.set_xticklabels(['1', '5,000', '10,000', '15,000', f'{max_rank:,}'])
            ax.set_yticks(np.linspace(0, 1, 5))
            if len(ax.lines) != 6:
                raise ValueError('Expected six CDF curves per panel')
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', ncol=3,
                   bbox_to_anchor=(0.5, 1), frameon=False, fontsize=10)
        fig.supxlabel('Rank among covered protein-coding genes (strongest → weakest)',
                      fontsize=12, y=0.02)
        fig.supylabel('Cumulative fraction of human genes', fontsize=12, x=0.004)
        fig.tight_layout(rect=(0.015, 0.06, 1, 0.9))
        for extension in ('png', 'pdf'):
            fig.savefig(OUT / f'human_gene_rank_panels_all_traits.{extension}', dpi=200)
        plt.close(fig)

    if sha256(notebook_path) != notebook_hash:
        raise ValueError('The analysis notebook changed during figure generation')
    report = {
        'gene_set': 'All 48 rows with Species == Human in ST6, header row 2',
        'symbol_mapping': 'Exact symbols; no aliases, orthology filters, or exclusions',
        'method': ('Same genome-wide covered-gene ranks as section 4E; minimum SNP p-value '
                   'in inclusive GENCODE TSS windows; strongest rank 1; average tied ranks; '
                   'empirical CDF across the 48 ST6 Human genes'),
        'validation': ('All saved genome-wide ranks independently checked by sorted-array '
                       'insertion positions; identical 48 human gene IDs in all 12 series'),
        'statistical_tests': 'None; this task creates descriptive CDFs only',
        'notebook_unchanged': True,
        'inputs_sha256': {str(p.relative_to(ROOT)): sha256(p)
                          for p in (workbook_path, ranks_path, notebook_path)},
        'environment': {'python': platform.python_version(), 'platform': platform.platform(),
                        'numpy': np.__version__, 'pandas': pd.__version__,
                        'matplotlib': matplotlib.__version__},
        'series': audit,
    }
    (OUT / 'human_gene_cdf_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(pd.DataFrame(audit).to_string(index=False))
    print('Validated all 48 Human genes in all 12 trait/window combinations.')
    print('Notebook SHA256 unchanged:', notebook_hash)
    print('Figure:', OUT / 'human_gene_rank_panels_all_traits.png')


if __name__ == '__main__':
    main()
