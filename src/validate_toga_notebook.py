"""Check TOGA notebook outputs, unchanged scoring inputs, and executed cells."""
import hashlib
import json
import platform
from pathlib import Path

import nbformat
import numpy as np
import pandas as pd
import scipy
from PIL import Image, ImageDraw, ImageOps
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'notebook/results/ipynb/toga'
notebook = nbformat.read(ROOT / 'dog_ocd_analysis.ipynb', as_version=4)
code = [c for c in notebook.cells if c.cell_type == 'code']
assert all(c.execution_count is not None for c in code)
assert not any(o.output_type == 'error' for c in code for o in c.outputs)
selected = pd.read_csv(OUT / 'toga_dog_gene_set.csv')
original = pd.read_excel(ROOT / 'data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx', sheet_name='ST6', header=1)
expected = original.loc[original.Species.eq('Dog') &
    original['TOGA (Human:Dog)**'].astype(str).str.strip().eq('1:1'), 'Gene*']
assert set(selected['Gene*']) == set(expected) and len(selected) == 83
ranked = pd.read_csv(OUT / 'mann_whitney_gene_ranks.csv')
tests = pd.read_csv(OUT / 'mann_whitney_results.csv')
assert len(tests) == 36
primary = tests.trait.eq('OCD')
assert primary.sum() == 6
np.testing.assert_allclose(tests.loc[primary, 'holm_adjusted_p'],
    multipletests(tests.loc[primary, 'one_sided_p'], method='holm')[1], rtol=1e-12, atol=1e-15)
assert tests.loc[~primary, 'holm_adjusted_p'].isna().all()
assert tests.loc[primary, 'holm_family_tests'].eq(6).all()
for row in tests.itertuples():
    frame = ranked.loc[ranked.trait.eq(row.trait) & ranked.window.eq(row.window)]
    np.testing.assert_array_equal(frame['rank'], frame.min_p.rank(method='average'))
    mask_col = {'combined': 'dog', 'significant': 'dog_significant',
                'suggestive': 'dog_suggestive'}[row.dog_set]
    chosen = frame[mask_col]
    assert chosen.sum() == {'combined': 83, 'significant': 43, 'suggestive': 40}[row.dog_set]
    assert set(frame.loc[frame.dog, 'gene_id']) == set(selected.human_gene_id)
    check = mannwhitneyu(frame.loc[chosen, 'min_p'], frame.loc[~frame.dog, 'min_p'],
                         alternative='less', method='asymptotic')
    np.testing.assert_allclose(check.pvalue, row.one_sided_p, rtol=1e-12, atol=1e-15)

for trait in tests.trait.unique():
    filename = f'per_gene_{trait.lower()}.csv'
    old_dir = ROOT / 'notebook/results' if trait in ['OCD', 'CAD'] else OUT.parent
    old = pd.read_csv(old_dir / filename, dtype={'chrom': str}).sort_values('gene_id').reset_index(drop=True)
    new = pd.read_csv(OUT / filename, dtype={'chrom': str}).sort_values('gene_id').reset_index(drop=True)
    columns = [c for c in old if not c.startswith('dog')]
    pd.testing.assert_frame_equal(old[columns], new[columns], check_dtype=False,
                                  check_exact=False, rtol=1e-12, atol=1e-15)
    assert set(new.loc[new.dog & ~old.dog, 'symbol']) == {
        'EMC2', 'PLXNA2', 'PPA2', 'PROSER3', 'SH3GL3'}
shifts = pd.read_csv(OUT / 'circular_shift_results.csv')
assert len(shifts) == 24 and shifts.Permutations.eq(10_000).all()
assert shifts.loc[shifts['Dog set'].eq('combined'), 'Dog genes'].eq(83).all()
assert shifts.loc[shifts['Dog set'].eq('significant'), 'Dog genes'].eq(43).all()
pd.testing.assert_frame_equal(pd.read_csv(OUT / 'human_mann_whitney_results.csv'),
                              pd.read_csv(OUT.parent / 'human_mann_whitney_results.csv'))
pd.testing.assert_frame_equal(pd.read_csv(OUT / 'human_gene_rank_cdf_source_data.csv'),
                              pd.read_csv(OUT.parent / 'human_gene_rank_cdf_source_data.csv'))

audit = {'orthology_source': 'Original ST6 TOGA 1:1 labels',
         'dog_genes': 83, 'significant': 43, 'suggestive': 40,
         'restored_genes': ['EMC2', 'PLXNA2', 'PPA2', 'PROSER3', 'SH3GL3'],
         'executed_code_cells': len(code), 'all_cells_error_free': True,
         'all_six_trait_scores_match_original_genome_wide_inputs': True,
         'mann_whitney_tests_independently_checked': 36,
         'dog_holm_family_tests': 6, 'dog_holm_traits': ['OCD'],
         'control_trait_holm_values_absent': True,
         'circular_shift_tests': 24, 'shifts_per_test': 10_000,
         'human_mann_whitney_and_cdf_data_unchanged': True,
         'environment': {'python': platform.python_version(), 'platform': platform.platform(),
                         'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
         'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in [ROOT / 'data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx',
                              ROOT / 'dog_ocd_analysis.ipynb', ROOT / 'src/toga_gene_set.py']}}
(OUT / 'toga_validation.json').write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps(audit, indent=2))
print(tests.loc[tests.dog_set.eq('combined'),
    ['trait', 'window', 'n_dog', 'one_sided_p', 'holm_adjusted_p']].to_string(index=False))

# Contact sheet for visual review; the scientific figures themselves are unchanged.
paths = sorted(OUT.glob('*.png')) + [ROOT / 'notebook/results/shet/shet_vs_ocd_500kb.png']
tile_width, tile_height, cols = 650, 420, 3
sheet = Image.new('RGB', (cols * tile_width, ((len(paths) + cols - 1) // cols) * tile_height), 'white')
draw = ImageDraw.Draw(sheet)
for i, p in enumerate(paths):
    picture = Image.open(p).convert('RGB')
    picture.thumbnail((tile_width - 12, tile_height - 32))
    x, y = (i % cols) * tile_width, (i // cols) * tile_height
    draw.text((x + 6, y + 3), p.name, fill='black')
    sheet.paste(picture, (x + (tile_width - picture.width) // 2, y + 25))
sheet.save('/private/tmp/dog_ocd_toga_figures.png')
