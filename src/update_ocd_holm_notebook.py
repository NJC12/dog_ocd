"""Update notebook correction policy, preserving the previous 36-test outputs."""
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
path = ROOT / 'dog_ocd_analysis.ipynb'
out = ROOT / 'notebook/results/ipynb/toga'
archive = out / 'holm36'
if not archive.exists():
    archive.mkdir()
    shutil.copy2(path, archive / 'executed_notebook.json')
    for p in out.iterdir():
        if p.is_file() and (p.name in ['METHODS.md', 'toga_validation.json',
                                      'psychiatric_independent_validation.json'] or
                            ('mann_whitney' in p.name and 'gene_ranks' not in p.name)):
            shutil.copy2(p, archive / p.name)

nb = json.loads(path.read_text())
def put(i, s):
    nb['cells'][i]['source'] = s.splitlines(keepends=True)

put(28, '''### Rank all covered genes

Calculate all 36 dog-set Mann–Whitney tests. The primary Holm family contains
only the **six OCD tests**: two window sizes × three dog sets. CAD, BIP, ADHD,
MDD and SCZ are controls and retain raw p-values; their Holm values are blank.
Repeated OCD tables show the same six tests. The CDF uses genome-wide ranks,
separately from the chromosome/SNP-count-stratified circular-shift ranks.
''')
put(32, '''#### Mann–Whitney tables

Each test compares a dog-gene set with other covered protein-coding genes
within the same trait and window. `P(dog stronger)` is `U / (n_dog × n_non-dog)`,
with ties contributing one half. Values above 0.5 indicate stronger dog-set signal.
Holm p-values apply only to the **six OCD tests** (two windows × three dog sets).
The five control traits report raw p-values and show a dash for Holm p.
Repeated OCD tables do not add tests to the correction family.
''')
put(48, '''#### Interpretation

The dog-set Holm family contains the six OCD tests. The five control traits
are reported with unadjusted p-values. Circular-shift p-values remain unadjusted.
Mann–Whitney tests do not account for overlapping windows, local LD, SNP-count
differences or clustering; circular shifts account for chromosome/SNP-count
strata and gene-order patterns. Different significance between traits does
not establish a difference in enrichment. The shared human OCD input in ST6
selection makes these findings exploratory.
''')
put(63, '''### D. Summary across all six traits

The table shows the combined 83-gene set at both windows. Holm correction is
restricted to the **six OCD tests**: two windows × three dog sets. The five
control traits retain raw p-values and have no reported Holm value. Repeated
OCD tables refer to the same six tests. Circular-shift p-values are unadjusted.
''')
put(67, '''#### Mann–Whitney tests

Each window table contains all three dog sets for six traits. Only OCD rows
have Holm p-values, from the six-test OCD family. Control traits show raw
p-values and a dash in the Holm column.
''')

for i, cell in enumerate(nb['cells']):
    s = ''.join(cell['source'])
    if cell['cell_type'] == 'code':
        s = s.replace('read_toga_gene_set, apply_dog_labels, dog_circular_shifts',
                      'read_toga_gene_set, apply_dog_labels, dog_circular_shifts, adjust_ocd_tests')
        s = s.replace("mw_results['holm_adjusted_p'] = multipletests(mw_results.one_sided_p, method='holm')[1]",
                      'mw_results = adjust_ocd_tests(mw_results)')
        s = s.replace('.style.format(formats)', ".style.format(formats, na_rep='—')")
        s = s.replace("'Non-dog genes': '{:,}'})", "'Non-dog genes': '{:,}'}, na_rep='—')")
        s = s.replace("'circular_shift_p': '{:.3g}'}).hide", "'circular_shift_p': '{:.3g}'}, na_rep='—').hide")
        if i == 64:
            s = s.replace('for trait in trait_frames:\n', '''for trait in trait_frames:
    if trait != 'OCD':
        display(Markdown(f"**{trait}:** Control trait; raw p-values reported without Holm correction."))
        continue
''')
        if i == 33:
            s = s.replace("f'{trait}: dog genes versus non-dog genes'",
                "f'{trait}: dog genes versus non-dog genes · Holm applies only to six OCD tests'")
        if i == 68:
            s = s.replace("f'Mann–Whitney tests · {WINDOW_LABELS[window]} TSS window'",
                "f'Mann–Whitney tests · {WINDOW_LABELS[window]} TSS window · Holm: six OCD tests only'")
        cell['outputs'] = []
        cell['execution_count'] = None
    else:
        s = s.replace('Holm adjustment spans all\n36 unique tests across the six traits.',
                      'Holm correction covers only the six OCD tests; this control trait reports raw p-values.')
        s = s.replace('Holm adjustment spans\n36 unique tests across all six traits.',
                      'Holm correction covers only the six OCD tests; this control trait reports raw p-values.')
        s = s.replace('separately from the 36 dog-set tests', 'separately from the six OCD dog-set tests')
    put(i, s)
path.write_text(json.dumps(nb, indent=1) + '\n')
print('Notebook updated; previous 36-test correction preserved in', archive)
